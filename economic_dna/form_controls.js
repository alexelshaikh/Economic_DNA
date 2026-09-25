(() => {
  const win = window;
  const doc = win.document;
  const previous = win.__dnaFormControls;
  if (previous) previous.dispose();
  const keys = Object.keys(config.committed);
  const fieldSelector = keys.map(key => `.st-key-${key}`).join(",");
  const calculateKeys = ["calculate_header", "calculate_scenario", "calculate_panel"];
  const root = key => doc.querySelector(`.st-key-${key}`);
  const input = key => root(key)?.querySelector("input");
  const read = key => {
    const field = input(key);
    if (!field) return undefined;
    if (field.type === "radio") {
      const checked = root(key).querySelector("input:checked");
      return checked?.closest("label").querySelector("p")?.textContent.trim();
    }
    if (field.type === "checkbox") return field.checked;
    if (field.type === "number") return field.value === "" ? NaN : Number(field.value);
    return field.value;
  };
  let frame = 0;
  let pending = 0;
  let actionValues = {};
  let queue = Promise.resolve();
  let disposed = false;
  const labelSidebarControls = () => {
    const open = doc.querySelector('[data-testid="stExpandSidebarButton"]');
    const close = doc.querySelector('[data-testid="stSidebarCollapseButton"] button');
    if (open) open.setAttribute("aria-label", "Open inputs");
    if (close) {
      close.setAttribute("aria-label", "Close inputs");
      const icon = close.querySelector('[data-testid="stIconMaterial"]');
      if (icon && icon.textContent !== "close") icon.textContent = "close";
    }
  };
  const refresh = () => {
    frame = 0;
    labelSidebarControls();
    let dirty = false;
    for (const key of keys) {
      const value = read(key);
      const changed = value !== undefined && value !== config.committed[key] && value !== config.displayed[key];
      root(key)?.classList.toggle("input-changed", changed);
      dirty ||= changed;
    }
    doc.documentElement.classList.toggle("inputs-pending", dirty);
    for (const key of calculateKeys) {
      const button = root(key)?.querySelector("button");
      const label = button?.querySelector("p");
      const text = dirty ? "Update charts" : "Calculate";
      if (label && label.textContent !== text) label.textContent = text;
      if (button) button.setAttribute("aria-label", dirty ? "Update charts. Inputs changed since the last calculation." : "Calculate");
    }
    const notice = doc.querySelector(".pending-notice");
    const message = dirty ? "Recalculate to update charts." : "Charts up to date";
    if (notice && notice.textContent !== message) notice.textContent = message;
  };
  const schedule = () => {
    if (!frame && !disposed) frame = win.requestAnimationFrame(refresh);
  };
  const settle = () => new Promise(resolve => win.requestAnimationFrame(() => win.setTimeout(resolve, 0)));
  const setValue = async (key, value) => {
    const field = input(key);
    if (read(key) === value) return;
    if (field.type === "radio") {
      const label = Array.from(root(key).querySelectorAll('[data-testid="stRadioOption"]'))
        .find(option => option.querySelector("p")?.textContent.trim() === value);
      label.click();
    } else if (field.type === "checkbox") {
      field.click();
    } else {
      // Use native input events so React and Streamlit's pending form state
      // both receive the edit, including for a currently hidden cost panel.
      Object.getOwnPropertyDescriptor(win.HTMLInputElement.prototype, "value").set.call(field, String(value));
      field.dispatchEvent(new win.Event("input", {bubbles: true}));
      await settle();
      field.dispatchEvent(new win.FocusEvent("focusout", {bubbles: true}));
    }
    await settle();
  };
  const actionFor = target => {
    if (!target.closest("button")) return null;
    return Object.entries(config.actions).find(([key]) => root(key)?.contains(target));
  };
  const click = event => {
    if (disposed) return;
    const submit = event.target.closest("button");
    if (pending && submit && [...calculateKeys, "theme_toggle"].some(key => root(key)?.contains(submit))) {
      event.preventDefault();
      event.stopImmediatePropagation();
      queue.then(() => submit.click());
      return;
    }
    const action = actionFor(event.target);
    if (!action) return;
    const [, operation] = action;
    const values = operation.values || {[operation.scale]: read(operation.scale) * operation.factor};
    if (Object.entries(values).some(([key, value]) => !input(key) || (typeof value === "number" && !Number.isFinite(value)))) return;
    event.preventDefault();
    event.stopImmediatePropagation();
    if (!pending) actionValues = {};
    pending += 1;
    queue = queue.then(async () => {
      try {
        const updates = operation.values || {
          [operation.scale]: (actionValues[operation.scale] ?? read(operation.scale)) * operation.factor,
        };
        Object.assign(actionValues, updates);
        await Promise.all(Object.entries(updates).map(([key, value]) => setValue(key, value)));
      } finally {
        pending -= 1;
        schedule();
      }
    });
  };
  const changed = event => {
    if (disposed) return;
    if (event.target.closest(fieldSelector)) schedule();
  };
  const keydown = event => {
    if (event.key !== "Escape" || !win.matchMedia("(max-width: 640px)").matches) return;
    const close = doc.querySelector('[data-testid="stSidebar"][aria-expanded="true"] [data-testid="stSidebarCollapseButton"] button');
    if (close) close.click();
  };
  // Delegation survives Streamlit reconciling widgets after Calculate or a
  // theme change. Only one observer/listener set is retained per document.
  doc.addEventListener("click", click, true);
  doc.addEventListener("input", changed, true);
  doc.addEventListener("change", changed, true);
  doc.addEventListener("keydown", keydown);
  const observer = new win.MutationObserver(records => {
    if (records.some(record => record.type === "childList" || (
      record.target.closest(fieldSelector)
      && record.oldValue !== record.target.getAttribute(record.attributeName)
    ))) schedule();
  });
  observer.observe(doc.body, {
    childList: true, subtree: true, attributes: true,
    attributeFilter: ["class", "value", "checked", "data-selected"], attributeOldValue: true,
  });
  win.__dnaFormControls = {
    dispose() {
      disposed = true;
      observer.disconnect();
      win.cancelAnimationFrame(frame);
      doc.removeEventListener("click", click, true);
      doc.removeEventListener("input", changed, true);
      doc.removeEventListener("change", changed, true);
      doc.removeEventListener("keydown", keydown);
    },
  };
  schedule();
})();
