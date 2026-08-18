#!/usr/bin/env python3
"""Plot the optimization history from benchmarks/history.json.

2x2 grid (delaunay/alpha x unweighted/weighted): one line per point-cloud
size, thick line = mean total_ms across the group's scenario variants
(two/distinct colors x half/quarter scale), shaded band = min..max across
those variants. Log y. Below: the latest profiled-function table from
profile_snapshot.json.

Usage: .venv/bin/python plot_history.py  ->  optimization_history.png
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).parent
BENCH = HERE.parent

GROUPS = {
    ("delaunay", "unweighted"): ["delaunay"],
    ("delaunay", "weighted"): ["delaunay_weighted_half", "delaunay_weighted_quarter"],
    ("alpha", "unweighted"): ["alpha_uniform_half", "alpha_uniform_quarter"],
    ("alpha", "weighted"): ["alpha_weighted_half", "alpha_weighted_quarter"],
}
SIZES = [1000, 5000, 20000, 50000]
SIZE_COLORS = {1000: "#9ecae1", 5000: "#6baed6", 20000: "#3182bd", 50000: "#08519c"}


def main() -> None:
    history = json.loads((BENCH / "history.json").read_text())
    labels = [e["commit_id"][:7] for e in history]
    x = range(len(history))

    fig = plt.figure(figsize=(15, 13))
    gs = GridSpec(3, 2, height_ratios=[3, 3, 2.1], hspace=0.45, wspace=0.18)

    for ax_i, ((pipeline, weighting), scenarios) in enumerate(GROUPS.items()):
        ax = fig.add_subplot(gs[ax_i // 2, ax_i % 2])
        for n in SIZES:
            series = []  # one list of totals per (scenario, colors) variant
            for scen in scenarios:
                for colors in ("two", "distinct"):
                    key = f"{n}/{colors}/{scen}"
                    vals = [e["tests"][key][1] for e in history if key in e["tests"]]
                    if len(vals) == len(history):
                        series.append(vals)
            if not series:
                continue
            mean = [sum(v[i] for v in series) / len(series) for i in x]
            lo = [min(v[i] for v in series) for i in x]
            hi = [max(v[i] for v in series) for i in x]
            c = SIZE_COLORS[n]
            ax.plot(x, mean, color=c, lw=2.8, label=f"n={n:,}", zorder=3)
            ax.fill_between(x, lo, hi, color=c, alpha=0.25, lw=0, zorder=2)
        ax.set_yscale("log")
        ax.set_title(f"{pipeline} / {weighting}", fontsize=12, fontweight="bold")
        ax.set_ylabel("total_ms (log)")
        ax.set_xticks(list(x))
        ax.set_xticklabels(labels, fontsize=6.5, rotation=60, ha="right")
        ax.set_xlabel("recorded iteration (commit)", fontsize=9)
        ax.grid(True, which="both", alpha=0.25)
        if ax_i == 0:
            ax.legend(loc="lower left", fontsize=9)

    # Profiled-function table
    snap = json.loads((HERE / "profile_snapshot.json").read_text())
    tax = fig.add_subplot(gs[2, :])
    tax.axis("off")
    cells = [[r["group"], f"{r['cpu_ms']}", f"{r['pct']:.1f}%"] for r in snap["rows"]]
    table = tax.table(
        cellText=cells,
        colLabels=["call group (self time)", "CPU-ms", "% of CPU"],
        colWidths=[0.78, 0.10, 0.12],
        cellLoc="left",
        loc="upper center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.0, 1.35)
    tax.set_title("Profiled main calls — " + snap["provenance"], fontsize=7.5, loc="left", pad=14, wrap=True)
    tax.text(0, -0.28, snap["note"], fontsize=7.5, style="italic", transform=tax.transAxes)

    fig.suptitle(
        "DelaunayInterfaces optimization history — interface-generation total_ms per iteration\n"
        "thick line: mean over scenario variants (two/distinct colors, half/quarter scale); band: min..max",
        fontsize=13,
    )
    out = HERE / "optimization_history.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
