#!/usr/bin/env python3
"""Genome synteny circos plot from a PAF file (pure matplotlib, no extra deps).

Reference chromosomes are laid over the top semicircle (left to right),
query chromosomes over the bottom semicircle, both proportional to length.
Syntenic blocks become filled bezier ribbons connecting the two rings,
colored by query chromosome; minus-strand blocks cross over (inversions).

Usage:
    python plot_syntenic_circos.py input.paf -o circos.png
    python plot_syntenic_circos.py input.paf -o circos.pdf --top 8 --min-anchors 50
"""
import argparse
import re

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Wedge

PALETTE = ["#4477AA", "#EE7733", "#228833", "#CCBB44", "#66CCEE",
           "#AA3377", "#AA4499", "#88CCEE", "#44AA99", "#999933"]
IDEO_REF_C = "#C7C7C7"

R_OUT, R_IN = 1.00, 0.93   # ideogram ring radii
R_LINK = 0.90              # ribbon anchor radius
R_CTRL = 0.15              # bezier control radius (pulls ribbons through centre)
R_LAB = 1.08               # label radius


def natural_key(name: str):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", name)]


def read_paf(path: str, min_anchors: int):
    links, ref_len, qry_len = [], {}, {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or line.startswith("#"):
                continue
            try:
                qs, qe, rs, re_ = int(f[2]), int(f[3]), int(f[7]), int(f[8])
                na = int(f[9]) if len(f) > 9 and f[9].lstrip("-").isdigit() else 1
            except ValueError:
                continue
            if na < min_anchors:
                continue
            links.append((f[0], qs, qe, f[4], f[5], rs, re_, na))
            qry_len[f[0]] = int(f[1])
            ref_len[f[5]] = int(f[6])
    return links, ref_len, qry_len


def pol(theta_deg: float, r: float):
    t = np.deg2rad(theta_deg)
    return r * np.cos(t), r * np.sin(t)


def cumulative(order, lengths):
    """cumulative bp offset per sequence"""
    cum, out = 0, {}
    for n in order:
        out[n] = cum
        cum += lengths[n]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paf")
    ap.add_argument("-o", "--out", default="syntenic_circos.png")
    ap.add_argument("--min-anchors", type=int, default=25)
    ap.add_argument("--top", type=int, default=0,
                    help="keep only the N longest query sequences (0 = all)")
    ap.add_argument("--max-links", type=int, default=3000,
                    help="draw at most this many largest ribbons")
    ap.add_argument("--alpha", type=float, default=0.55)
    ap.add_argument("--title", default="")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    links, ref_len, qry_len = read_paf(args.paf, args.min_anchors)
    if not links:
        raise SystemExit("no alignment blocks in input")

    ref_order = sorted(ref_len, key=natural_key)
    qry_order = sorted(qry_len, key=natural_key)
    if args.top and len(qry_order) > args.top:
        keep = set(sorted(qry_order, key=lambda n: -qry_len[n])[: args.top])
        qry_order = [n for n in qry_order if n in keep]
        links = [l for l in links if l[0] in keep]

    ref_cum = cumulative(ref_order, ref_len)
    qry_cum = cumulative(qry_order, qry_len)
    ref_dpb = 180.0 / sum(ref_len[n] for n in ref_order)
    qry_dpb = 180.0 / sum(qry_len[n] for n in qry_order)
    qry_color = {n: PALETTE[i % len(PALETTE)] for i, n in enumerate(qry_order)}

    # angle maps: ref covers the TOP semicircle 180 -> 0 degrees (left to
    # right), query the BOTTOM semicircle 180 -> 360
    def ang_ref(name, pos):
        return 180.0 - (ref_cum[name] + pos) * ref_dpb

    def ang_qry(name, pos):
        return 180.0 + (qry_cum[name] + pos) * qry_dpb

    # keep only the largest links for drawing clarity
    links.sort(key=lambda l: -(min(l[6] - l[5], l[2] - l[1])))
    links = links[: args.max_links]

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect("equal")
    ax.axis("off")
    # patches do not set data limits: fix the canvas so the circle fills
    # the square figure
    lim = 1.30
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    def draw_ideogram(name, ang0, ang1, length, color):
        lo, hi = min(ang0, ang1), max(ang0, ang1)
        ax.add_patch(Wedge((0, 0), R_OUT, lo, hi, width=R_OUT - R_IN,
                           facecolor=color, edgecolor="white", lw=1.0, zorder=3))
        # Mb ticks (no labels), at most ~6 per chromosome
        step = max(1e6, round(length / 6 / 1e6) * 1e6)
        p = step
        while p < length:
            tt = ang0 + (ang1 - ang0) * (p / length)
            x0, y0 = pol(tt, R_OUT)
            x1, y1 = pol(tt, R_OUT + 0.018)
            ax.plot([x0, x1], [y0, y1], color="#909090", lw=0.6, zorder=4)
            p += step
        # radial label
        mid = (ang0 + ang1) / 2
        m = mid % 360
        rot = m - 180 if 90 < m < 270 else m
        ax.text(*pol(mid, R_LAB), name, rotation=rot, rotation_mode="anchor",
                ha="center", va="center", fontsize=8.5, color="#333333", zorder=5)

    for n in ref_order:
        draw_ideogram(n, ang_ref(n, 0), ang_ref(n, ref_len[n]), ref_len[n], IDEO_REF_C)
    for n in qry_order:
        draw_ideogram(n, ang_qry(n, 0), ang_qry(n, qry_len[n]), qry_len[n], qry_color[n])

    for qn, qs, qe, strand, rn, rs, re_, _na in links:
        a1, a2 = ang_ref(rn, rs), ang_ref(rn, re_)
        b1, b2 = ang_qry(qn, qs), ang_qry(qn, qe)
        # minus-strand blocks swap endpoints -> the ribbon crosses over
        if strand == "-":
            b1, b2 = b2, b1
        A1, A2 = pol(a1, R_LINK), pol(a2, R_LINK)
        B1, B2 = pol(b1, R_LINK), pol(b2, R_LINK)
        cA1, cA2 = pol(a1, R_CTRL), pol(a2, R_CTRL)
        cB1, cB2 = pol(b1, R_CTRL), pol(b2, R_CTRL)
        verts = [A1, cA1, cB1, B1, B2, cB2, cA2, A2, A1]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                 Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                 Path.CLOSEPOLY]
        ax.add_patch(PathPatch(Path(verts, codes), facecolor=qry_color[qn],
                               edgecolor="none", alpha=args.alpha, zorder=2))

    ax.text(0, 0.03, args.title or "Genome synteny", ha="center", va="center",
            fontsize=14, color="#222222", zorder=5)
    ax.text(0, -0.05, f"{len(links)} syntenic blocks | ribbons colored by query chromosome",
            ha="center", va="center", fontsize=8.5, color="#666666", zorder=5)

    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    print(f"wrote {args.out}  ({len(links)} ribbons, {len(ref_order)} ref / {len(qry_order)} query sequences)")


if __name__ == "__main__":
    main()
