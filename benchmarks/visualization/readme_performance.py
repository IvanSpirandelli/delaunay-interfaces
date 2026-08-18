#!/usr/bin/env python3
"""README performance figure: current speed of interface generation.

Left: time to compute one interface vs point-cloud size on FOUR-COLOR random
clouds (thick line = mean over the measured runs, band = min..max over runs).
Data comes from bench_results.json, produced by readme_bench.cpp.
Right: where one computation spends its CPU time, per pipeline (sampled
profile of compute_interface_surface only - no benchmark harness).

Usage: .venv/bin/python readme_performance.py -> readme_performance.png
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).parent

GROUP_COLORS = {
    "full Delaunay, unweighted": "#08519c",
    "full Delaunay, weighted": "#3182bd",
    "alpha complex, weighted": "#e6550d",
    "alpha complex, uniform radius": "#fd8d3c",
}
# One fixed color per profile category, identical across both pies.
CATEGORY_COLORS = ["#6baed6", "#e6550d", "#3182bd", "#08519c", "#fd8d3c", "#cccccc"]

plt.rcParams.update({
    "font.size": 13,
    "axes.titlesize": 15,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
})


def main() -> None:
    bench = json.loads((HERE / "bench_results.json").read_text())
    reps = bench["reps"]

    fig = plt.figure(figsize=(15, 7))
    gs = GridSpec(2, 2, width_ratios=[1.6, 1.0], hspace=0.55, wspace=0.3)

    ax = fig.add_subplot(gs[:, 0])
    for label, runs_by_n in bench["groups"].items():
        sizes = sorted(int(n) for n in runs_by_n)
        mean = [sum(runs_by_n[str(n)]) / len(runs_by_n[str(n)]) for n in sizes]
        lo = [min(runs_by_n[str(n)]) for n in sizes]
        hi = [max(runs_by_n[str(n)]) for n in sizes]
        c = GROUP_COLORS[label]
        ax.plot(sizes, mean, color=c, lw=3.0, marker="o", ms=6, label=label, zorder=3)
        ax.fill_between(sizes, lo, hi, color=c, alpha=0.25, lw=0, zorder=2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(sizes)
    ax.set_xticklabels([f"{n:,}" for n in sizes])
    ax.minorticks_off()
    ax.set_xlabel("points in the cloud")
    ax.set_ylabel("time to compute the interface (ms)")
    ax.set_title("Interface generation time", fontweight="bold")
    ax.grid(True, which="major", alpha=0.3)
    ax.legend(loc="upper left")
    note = f"mean of {reps} runs, band: min..max over runs\nApple M-series laptop"
    ax.text(0.98, 0.02, note, transform=ax.transAxes, fontsize=10, ha="right", va="bottom", alpha=0.85)

    snap = json.loads((HERE / "profile_snapshot.json").read_text())
    categories = snap["categories"]
    for i, (title, vals) in enumerate(snap["pipelines"].items()):
        pax = fig.add_subplot(gs[i, 1])
        pax.pie(
            vals,
            labels=None,
            colors=CATEGORY_COLORS,
            startangle=90,
            counterclock=False,
            autopct=lambda p: f"{p:.0f}%" if p >= 5 else "",
            textprops={"fontsize": 11},
            wedgeprops={"lw": 0.8, "ec": "white"},
            radius=1.15,
        )
        if i == 0:
            pax.legend(categories, fontsize=10, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)
        pax.set_title(f"Where the time goes — {title}\n({snap['input']})", fontsize=12)

    fig.suptitle("delaunay-interfaces: four colored point cloud in, interface surface out", fontsize=17, y=1.01)
    fig.text(0.985, 0.0, snap["note"], fontsize=8.5, style="italic", ha="right")
    out = HERE / "readme_performance.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
