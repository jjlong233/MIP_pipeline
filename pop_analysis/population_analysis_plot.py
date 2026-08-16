#!/usr/bin/env python3
"""
population_analysis_plot.py 

This script is for generating following three outputs:
  (A) PCA scatter (PC1 vs PC2) with 95% confidence ellipses per population
  (B) ADMIXTURE barplots for K = 1..4, individuals grouped by population
  (C) cross-validation error vs K  (optional extra panel, --cv-panel)

Usage:
    python population_analysis_plot.py --dir popstructure --popmap pop_map.txt --out Fig6_new

Inputs expected in --dir (produced by run_popstructure.sh):
    mip.pca.eigenvec, mip.pca.eigenval, mip.admix.<K>.Q, mip.admix.fam,
    cv_error.tsv
"""
import argparse
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# population order + colours (HP, LR, RU, YR as in the manuscript legend)
POP_ORDER = ["YR", "LR", "RU", "HP"]
POP_COLOR = {"YR": "#2C7FB8", "LR": "#41B6C4", "RU": "#F0A202", "HP": "#D7301F"}
POP_LABEL = {"YR": "YR (Yellow River)", "LR": "LR (Liao River)",
             "RU": "RU (Vladivostok)", "HP": "HP (E. hepuensis)"}


def read_popmap(path):
    pop = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                pop[parts[0]] = parts[1]
    return pop


def read_eigen(d):
    """PLINK .eigenvec (FID IID PC1..PCn) and .eigenval."""
    ids, pcs = [], []
    with open(os.path.join(d, "mip.pca.eigenvec")) as fh:
        for line in fh:
            p = line.split()
            if not p:
                continue
            ids.append(p[1])                      # IID
            pcs.append([float(x) for x in p[2:]])
    ev = []
    with open(os.path.join(d, "mip.pca.eigenval")) as fh:
        for line in fh:
            line = line.strip()
            if line:
                ev.append(float(line))
    return ids, np.array(pcs), np.array(ev)


def confidence_ellipse(x, y, ax, n_std=2.4477, **kw):
    """
    95% confidence ellipse for a bivariate normal (matches ggplot2
    stat_ellipse defaults: level=0.95, multivariate t/normal radius).
    n_std = sqrt(qchisq(0.95, df=2)) = 2.4477 for the normal approximation.
    """
    if x.size < 3:
        return None
    cov = np.cov(x, y)
    if not np.all(np.isfinite(cov)) or np.linalg.det(cov) <= 0:
        return None
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    theta = np.degrees(np.arctan2(*vecs[:, 0][::-1]))
    w, h = 2 * n_std * np.sqrt(vals)
    e = Ellipse((x.mean(), y.mean()), width=w, height=h, angle=theta, **kw)
    ax.add_patch(e)
    return e


def panel_pca(ax, ids, pcs, ev, pop):
    var = 100.0 * ev / ev.sum()
    for p in POP_ORDER:
        idx = [i for i, s in enumerate(ids) if pop.get(s) == p]
        if not idx:
            continue
        x, y = pcs[idx, 0], pcs[idx, 1]
        c = POP_COLOR[p]
        ax.scatter(x, y, s=42, color=c, edgecolor="black", linewidth=0.4,
                   label=POP_LABEL[p], zorder=3)
        confidence_ellipse(x, y, ax, facecolor=c, alpha=0.13,
                           edgecolor=c, linewidth=1.2, zorder=1)
    ax.axhline(0, color="grey", lw=0.5, ls=":", zorder=0)
    ax.axvline(0, color="grey", lw=0.5, ls=":", zorder=0)
    ax.set_xlabel(f"PC1 ({var[0]:.2f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.2f}%)")
    ax.set_title("A", loc="left", fontweight="bold", fontsize=13)
    ax.legend(frameon=False, fontsize=8, loc="best")
    return var


def panel_admixture(axes, d, order, pop, ks=(1, 2, 3, 4)):
    # individuals grouped by population, stable order within group
    idx = sorted(range(len(order)),
                 key=lambda i: (POP_ORDER.index(pop.get(order[i], "YR")), order[i]))
    labels = [order[i] for i in idx]
    pops = [pop.get(s, "?") for s in labels]

    for ax, K in zip(axes, ks):
        qf = os.path.join(d, f"mip.admix.{K}.Q")
        if not os.path.exists(qf):
            ax.set_visible(False)
            continue
        Q = np.loadtxt(qf, ndmin=2)[idx, :]
        bottom = np.zeros(Q.shape[0])
        palette = plt.cm.Set2(np.linspace(0, 1, max(Q.shape[1], 3)))
        for k in range(Q.shape[1]):
            ax.bar(range(Q.shape[0]), Q[:, k], bottom=bottom, width=1.0,
                   color=palette[k], edgecolor="none")
            bottom += Q[:, k]
        ax.set_xlim(-0.5, Q.shape[0] - 0.5)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 0.5, 1])
        ax.set_ylabel(f"K={K}", fontsize=9)
        ax.set_xticks([])
        for s in ("top", "right", "bottom"):
            ax.spines[s].set_visible(False)

    # population dividers + labels under the bottom panel
    last = axes[-1]
    bounds, start = [], 0
    for i in range(1, len(pops) + 1):
        if i == len(pops) or pops[i] != pops[start]:
            bounds.append((pops[start], start, i))
            start = i
    for name, s, e in bounds:
        last.text((s + e - 1) / 2.0, -0.22, name, ha="center", va="top",
                  fontsize=9, transform=last.get_xaxis_transform())
        if e < len(pops):
            for ax in axes:
                if ax.get_visible():
                    ax.axvline(e - 0.5, color="black", lw=0.8)
    axes[0].set_title("B", loc="left", fontweight="bold", fontsize=13)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="popstructure")
    ap.add_argument("--popmap", default="pop_map.txt")
    ap.add_argument("--out", default="Fig6_new")
    ap.add_argument("--ks", default="1,2,3,4")
    ap.add_argument("--cv-panel", action="store_true",
                    help="add a CV-error vs K panel")
    args = ap.parse_args()

    pop = read_popmap(args.popmap)
    ids, pcs, ev = read_eigen(args.dir)

    missing = [s for s in ids if s not in pop]
    if missing:
        sys.stderr.write(f"WARNING: {len(missing)} sample(s) absent from popmap: "
                         f"{', '.join(missing[:6])}\n")

    fam = os.path.join(args.dir, "mip.admix.fam")
    order = ([l.split()[1] for l in open(fam)] if os.path.exists(fam) else ids)
    ks = [int(k) for k in args.ks.split(",")]

    nrow = 1 + len(ks) + (1 if args.cv_panel else 0)
    fig = plt.figure(figsize=(8.2, 4.2 + 0.62 * len(ks) + (1.6 if args.cv_panel else 0)))
    gs = fig.add_gridspec(nrow, 1,
                          height_ratios=[4.2] + [0.62] * len(ks) +
                                        ([1.6] if args.cv_panel else []),
                          hspace=0.32)

    ax_pca = fig.add_subplot(gs[0, 0])
    var = panel_pca(ax_pca, ids, pcs, ev, pop)

    ax_adm = [fig.add_subplot(gs[1 + i, 0]) for i in range(len(ks))]
    panel_admixture(ax_adm, args.dir, order, pop, ks=ks)

    if args.cv_panel:
        axc = fig.add_subplot(gs[-1, 0])
        cvf = os.path.join(args.dir, "cv_error.tsv")
        if os.path.exists(cvf):
            arr = np.loadtxt(cvf, ndmin=2)
            axc.plot(arr[:, 0], arr[:, 1], "o-", color="black", ms=4)
            best = arr[np.argmin(arr[:, 1]), 0]
            axc.axvline(best, color="red", ls="--", lw=1)
            axc.set_xlabel("K"); axc.set_ylabel("CV error")
            axc.set_title("C", loc="left", fontweight="bold", fontsize=13)

    for ext in ("png", "pdf"):
        fig.savefig(f"{args.out}.{ext}", dpi=350, bbox_inches="tight")
    print(f"wrote {args.out}.png / {args.out}.pdf")
    print(f"PC1 {var[0]:.2f}% , PC2 {var[1]:.2f}% , PC3 {var[2]:.2f}%")


if __name__ == "__main__":
    main()
