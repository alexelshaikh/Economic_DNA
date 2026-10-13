"""Evaluate F1's displayed trajectories without refitting or changing observations."""

DIAGNOSTICS = '''
# Evaluate at observation dates, including fractional sequencing years, not future projections.
fit_observations = observations.copy()
fit_rows = []
for process, price, decline in (
    ("Synthesis", BASE.dna_synthesis_cost_per_mb, BASE.synthesis_decline_percent),
    ("Sequencing", BASE.dna_sequencing_cost_per_mb, BASE.sequencing_decline_percent),
):
    mask = fit_observations.process == process
    points = fit_observations.loc[mask]
    actual = points.cost_usd_per_logical_mb.to_numpy()
    predicted = price * (1 - decline / 100) ** (points.year.to_numpy() - BASE.dna_cost_base_year)
    residual = actual - predicted
    sse = float(np.sum(residual ** 2))
    sst = float(np.sum((actual - actual.mean()) ** 2))
    r_squared = 1 - sse / sst if sst > 0 else np.nan
    rmse = float(np.sqrt(np.mean(residual ** 2)))
    fit_observations.loc[mask, "predicted_usd_per_logical_mb"] = predicted
    fit_observations.loc[mask, "residual_usd_per_logical_mb"] = residual
    fit_rows.append({"process": process, "n_observations": len(points),
                     "first_year": float(points.year.min()), "last_year": float(points.year.max()),
                     "r_squared": r_squared, "rmse_usd_per_logical_mb": rmse})
f1_fit_metrics = pd.DataFrame(fit_rows)
print("F1 goodness of fit: displayed model trajectories, no refitting.")
print("R^2 and RMSE use historical observations on the original (non-logarithmic) cost scale.")
for row in f1_fit_metrics.itertuples():
    print(f"{row.process}: R^2 = {row.r_squared:.6f}; RMSE = {row.rmse_usd_per_logical_mb:,.6f} "
          f"USD/logical MB (n = {row.n_observations}, years {row.first_year:.2f}-{row.last_year:.2f}).")
f1_metrics_caption = " Goodness of fit against historical observations, on the original USD/logical MB scale: " + "; ".join(
    f"{row.process}: R2={row.r_squared:.4f}, RMSE={row.rmse_usd_per_logical_mb:.6g} USD/logical MB "
    f"(n={row.n_observations})" for row in f1_fit_metrics.itertuples()) + ". The plotted trajectories are not refitted."
if not SHOW_RECORDED_COSTS:
    observations = observations.iloc[:0].copy()
'''


def add_f1_diagnostics(source):
    if "f1_fit_metrics =" in source:
        return source
    start = source.index("if SHOW_RECORDED_COSTS:\n")
    plot_start = source.index('    for process, marker in', start)
    load = source[start + len("if SHOW_RECORDED_COSTS:\n"):plot_start]
    load = "\n".join(line[4:] if line.startswith("    ") else line for line in load.split("\n"))
    source = source[:start] + load + DIAGNOSTICS + "\nif SHOW_RECORDED_COSTS:\n" + source[plot_start:]
    source = source.replace('else "Historical points are disabled."),',
                            'else "Historical points are disabled.") + f1_metrics_caption,', 1)
    source = source.replace('"observations": observations}',
                            '"observations": observations, "fit_metrics": f1_fit_metrics, "fit_observations": fit_observations}', 1)
    return source
