import math

import numpy as np
from sympy import (
    Symbol,
    Function,
    Expr,
    summation,
    simplify,
    Piecewise,
    Add, Sum
)
from sympy.core.function import UndefinedFunction
from numbers import Number

class CostModel:
    _counter = 0

    def __init__(self, model_name=None):
        if model_name is None:
            CostModel._counter += 1
            model_name = f"Model_{CostModel._counter}"
        self.model_name = model_name

        # —––––– Core “static” symbols ––––––
        # n = number of objects (integer ≥ 0)
        self.n = Symbol(f"n_{model_name}", integer=True, nonnegative=True)

        # k = number of reads per period (integer ≥ 0)
        self.k = Symbol(f"k_{model_name}", integer=True, nonnegative=True)

        # obj_size = size per object (MB, >0)
        self.obj_size = Symbol(f"obj_size_{model_name}", positive=True)

        # device_durability(t): by default a Symbol; user can override with a Python function
        #self.device_durability = Symbol(f"device_durability_{model_name}")
        self.device_durability = Function(f"device_durability_{model_name}")

        # start_year = first year of operation (integer)
        self.start_year = Symbol(f"start_year_{model_name}", integer=True)

        # d = total duration (in years) of the analysis (integer ≥ 0)
        self.d = Symbol(f"d_{model_name}", integer=True, nonnegative=True)

        # max_repls = maximum number of replacements to include (integer ≥ 0)
        self.max_repls = Symbol(f"max_replacements_count_{model_name}",
                                integer=True, nonnegative=True)

        # —––––– cost‐per‐MB placeholders as Sympy Functions of t ––––––
        # Default each one to be a Function(t).  If you want time‐constant cost, override
        # with a Symbol or a Sympy Expr in t.
        self.cost_write = Function(f"cost_write_{model_name}")   # cost_write(t)
        self.cost_read  = Function(f"cost_read_{model_name}")    # cost_read(t)
        self.cost_maint = Function(f"cost_maint_{model_name}")   # cost_maint(t)

        # —––––– Y(i): an abstract “i-th replacement year” function ––––––
        # We will only need this if device_durability is non‐callable or symbolic; but
        # if device_durability is Python‐callable, we compute Y(i) numerically and do not
        # leave it abstract.  Still, we define Y(i) so that the expression is well‐formed.
        self.Y = Function(f"Y_{model_name}")

    def symbols_dict(self):
        """
        Return all the core Symbols/Functions in this model.
        """
        return {
            "n": self.n,
            "k": self.k,
            "obj_size": self.obj_size,
            "device_durability": self.device_durability,
            "start_year": self.start_year,
            "d": self.d,
            "max_repls": self.max_repls,
            "cost_write":   self.cost_write,
            "cost_read":    self.cost_read,
            "cost_maint":   self.cost_maint,
            "Y": self.Y
        }

    def _get_initial_write_cost(self):
        """
        initial_write =
          – if cost_write is a Sympy Function(t), cost_write(start_year)*n*obj_size
          – if cost_write is a Sympy Expr in t, substitute t→start_year
          – otherwise (Symbol or number), cost_write * n * obj_size
        """
        cw = self.cost_write

        t  = Symbol("t", integer=True)

        if isinstance(cw, UndefinedFunction):
            return cw(self.start_year) * self.n * self.obj_size

        if isinstance(cw, Expr) and (t in cw.free_symbols):
            return cw.subs(t, self.start_year) * self.n * self.obj_size

        # time‐constant Symbol or numeric
        return cw * self.n * self.obj_size

    def _get_maintenance_sum(self):
        """
        maintenance_sum =
          ∑_{t = start_year}^{start_year + d} [ cost_maint(t) * n * obj_size ].
        """
        cm = self.cost_maint
        t = Symbol("t", integer=True)

        # 1) If cost_maint is an “undefined Sympy Function” (e.g. Function("f")), call it as cm(t):
        if isinstance(cm, UndefinedFunction):
            term = cm(t)

        # 2) Otherwise (whether cm is a Sympy Expr _in_ t, or a Symbol/number _not_ in t),
        #    we can just treat cm itself as the “summand.”  If cm depends on t, Sympy will
        #    sum it correctly.  If cm does _not_ depend on t, Sympy will multiply by (d+1) for us.
        else:
            term = cm

        return summation(term * self.n * self.obj_size,
                         (t, self.start_year, self.start_year + self.d))

    def _get_read_sum(self):
        """
        read_sum =
          ∑_{t = start_year}^{start_year + d} [ cost_read(t) * k * obj_size ].
        """
        cr = self.cost_read
        t = Symbol("t", integer=True)

        # 1) If cost_read is an undefined Sympy Function, call cr(t)
        if isinstance(cr, UndefinedFunction):
            term = cr(t)
        else:
            # 2) Otherwise, cr can be:
            #    – a Sympy Expr already containing t, or
            #    – a Symbol/number (time‐constant).
            term = cr

        return summation(
            term * self.k * self.obj_size,
            (t, self.start_year, self.start_year + self.d)
        )

    def _get_replacement_cost(self):
        """
        Optimized replacement-cost calculation: unified numeric vs symbolic handling.
        Supports device_durability as Python callable, Sympy expression, or constant.
        """
        rd = self.device_durability
        cw = self.cost_write
        t = Symbol('t', integer=True)

        def cost_term(year):
            # Build cost_write(year) * n * obj_size
            if isinstance(cw, UndefinedFunction) or callable(cw):
                return cw(year) * self.n * self.obj_size
            elif isinstance(cw, Expr) and t in cw.free_symbols:
                return cw.subs(t, year) * self.n * self.obj_size
            else:
                return cw * self.n * self.obj_size

        # Attempt to get an integer for max_repls
        try:
            max_repls = int(self.max_repls)
        except Exception:
            max_repls = None

        # Determine if start_year or d are symbolic
        start_sym = isinstance(self.start_year, Expr) and bool(self.start_year.free_symbols)
        d_sym = isinstance(self.d, Expr) and bool(self.d.free_symbols)
        numeric_window = (max_repls is not None and not start_sym and not d_sym)

        # Treat Python numbers as valid constant durability
        is_rd_constant = isinstance(rd, Number) or (isinstance(rd, Expr) and not rd.free_symbols)
        is_rd_valid = callable(rd) or isinstance(rd, Expr) or isinstance(rd, Number)

        # Numeric branch: build actual replacement years list if rd yields numeric lifetimes
        if numeric_window and is_rd_valid:
            start_val = float(self.start_year)
            duration = float(self.d)
            window_end = start_val + duration
            years = []
            prev = start_val
            # determine constant lifetime if rd is constant
            const_life = None
            if is_rd_constant:
                const_life = float(rd) if isinstance(rd, Number) else float(rd)
            for _ in range(max_repls):
                # get lifetime
                if const_life is not None:
                    lifetime = const_life
                else:
                    lifetime_raw = rd(prev) if callable(rd) else rd.subs(t, prev)
                    try:
                        lifetime = float(lifetime_raw)
                    except Exception:
                        # cannot convert lifetime; bail to symbolic
                        years = []
                        break
                nxt = prev + lifetime
                if nxt > window_end:
                    break
                years.append(nxt)
                prev = nxt
            if years:
                return Add(*[cost_term(y) for y in years])
            # if no valid years, fall through to symbolic branch

        # Symbolic branch: build piecewise terms per potential replacement
        if max_repls is not None and is_rd_valid:
            pieces = []
            prev = self.start_year
            for _ in range(max_repls):
                if is_rd_constant:
                    lifetime_raw = rd
                else:
                    lifetime_raw = rd(prev) if callable(rd) else rd.subs(t, prev)
                next_expr = simplify(prev + lifetime_raw)
                cond = next_expr <= (self.start_year + self.d)
                pieces.append(Piecewise((cost_term(next_expr), cond), (0, True)))
                prev = next_expr
            return Add(*pieces)

        # Fallback placeholder if building schedule not possible
        return Symbol(f"replacement_costs_{self.model_name}")

    def eval(self):
        """
        Build the total cost:
          total_cost =
             initial_write
           + maintenance_sum
           + read_sum
           + replacement_cost

        Returns a symbolic expression.  If you want the “automatic”
        replacement‐years calculation, you must substitute numeric values for
        {start_year, d, max_repls} (and supply a Python function for
        device_durability) before evaluating.
        """
        initial_write    = self._get_initial_write_cost()
        maintenance_sum  = self._get_maintenance_sum()
        read_sum         = self._get_read_sum()
        replacement_cost = self._get_replacement_cost()

        total = initial_write + maintenance_sum + read_sum + replacement_cost
        return simplify(total)
