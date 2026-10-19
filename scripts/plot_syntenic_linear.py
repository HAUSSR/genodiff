#!/usr/bin/env python3
"""Genome synteny plot (linear, two stacked rows) from a PAF file.

Reference chromosomes form the top row, query chromosomes the bottom row,
both horizontal and proportional to length. Syntenic blocks become filled
bezier ribbons between the rows, colored by query chromosome; minus-strand
blocks cross over (inversions). Pure matplotlib.

Output .png/.pdf renders a plain figure; output .svg renders a hybrid SVG:
900-dpi PNG base with all text as separately editable, exactly anchored
<text> elements.

Usage:
    python plot_syntenic_linear.py input.paf -o synteny.png
    python plot_syntenic_linear.py input.paf -o synteny.svg --top 8 --min-anchors 50
"""
import argparse
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import PathPatch, Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svg_export

PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7",
           "#56B4E9", "#B07AA1", "#8C564B", "#999999", "#117733"]
IDEO_REF_FC = "#D2D2D2"
IDEO_REF_EC = "#7A7A7A"
NAME_C = "#1A1A1A"
MB_C = "#666666"

BAR_H = 0.085          # chromosome bar height (axes fraction)
Y_TOP_BAR = 0.735      # top row bar lower edge
Y_BOT_BAR = 0.205      # bottom row bar lower edge
Y_MID = (Y_TOP_BAR + Y_BOT_BAR + BAR_H) / 2  # ribbon bezier control height
X0, X1 = 0.02, 0.985   # row span (axes fraction)


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


def row_layout(order, lengths, equal):
    """fractional span per sequence + bp->fraction converter"""
    if equal:
        n = len(order)
        span = (X1 - X0) / n
        off = {name: X0 + i * span for i, name in enumerate(order)}
        def to_frac(name, pos):
            return off[name] + min(pos / max(lengths[name], 1), 1.0) * span
        return off, span, to_frac
    total = sum(lengths[n] for n in order)
    cum, off = 0.0, {}
    for n in order:
        off[n] = X0 + (X1 - X0) * cum / total
        cum += lengths[n]
    scale = (X1 - X0) / total
    def to_frac(name, pos):
        return off[name] + pos * scale
    return off, None, to_frac


def nice_step(length: float, max_ticks: int = 6):
    mb = 1e6
    for s in (0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100):
        if length / (s * mb) <= max_ticks:
            return s * mb
    return 200 * mb


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paf")
    ap.add_argument("-o", "--out", default="synteny_linear.png")
    ap.add_argument("--min-anchors", type=int, default=25)
    ap.add_argument("--top", type=int, default=0,
                    help="keep only the N longest query sequences (0 = all)")
    ap.add_argument("--max-links", type=int, default=3000,
                    help="draw at most this many largest ribbons")
    ap.add_argument("--equal", action="store_true",
                    help="equal width per chromosome (nice with many scaffolds)")
    ap.add_argument("--alpha", type=float, default=0.6)
    ap.add_argument("--title", default="")
    ap.add_argument("--dpi", type=int, default=300,
                    help="raster dpi for .png/.pdf output (SVG base is always 900)")
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

    ref_off, ref_span, ref_frac = row_layout(ref_order, ref_len, args.equal)
    qry_off, qry_span, qry_frac = row_layout(qry_order, qry_len, args.equal)
    qry_color = {n: PALETTE[i % len(PALETTE)] for i, n in enumerate(qry_order)}

    links.sort(key=lambda l: -(min(l[6] - l[5], l[2] - l[1])))
    links = links[: args.max_links]

    fig, ax = plt.subplots(figsize=(14, 6))
    ax.set_position([0.01, 0.01, 0.985, 0.98])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # reference row (top): gray bars, names above
    tr_x = ax.get_xaxis_transform()  # x in data (fraction), y in axes fraction
    for name in ref_order:
        x0 = ref_off[name]
        w = (ref_span if args.equal else (ref_frac(name, ref_len[name]) - x0))
        ax.add_patch(Rectangle((x0, Y_TOP_BAR), w, BAR_H, fc=IDEO_REF_FC,
                               ec=IDEO_REF_EC, lw=0.8, zorder=3,
                               transform=tr_x))
        # Mb ticks + numbers above the bar
        step = nice_step(ref_len[name])
        p = step
        while p < ref_len[name] * 0.999:
            x = ref_frac(name, p)
            ax.plot([x, x], [Y_TOP_BAR + BAR_H, Y_TOP_BAR + BAR_H + 0.012],
                    color=IDEO_REF_EC, lw=0.6, zorder=3, transform=tr_x)
            ax.text(x, Y_TOP_BAR + BAR_H + 0.018, f"{p / 1e6:g}", ha="center",
                    va="bottom", fontsize=6.5, color=MB_C,
                    transform=tr_x, rotation_mode="anchor")
            p += step
        ax.text(x0 + w / 2, Y_TOP_BAR + BAR_H + 0.052, name, ha="center",
                va="bottom", fontsize=9.5, color=NAME_C, fontweight="bold",
                transform=tr_x, rotation_mode="anchor")

    # query row (bottom): colored bars, names below
    for name in qry_order:
        x0 = qry_off[name]
        w = (qry_span if args.equal else (qry_frac(name, qry_len[name]) - x0))
        ax.add_patch(Rectangle((x0, Y_BOT_BAR), w, BAR_H, fc=qry_color[name],
                               ec="#333333", lw=0.6, zorder=3, transform=tr_x))
        step = nice_step(qry_len[name])
        p = step
        while p < qry_len[name] * 0.999:
            x = qry_frac(name, p)
            ax.plot([x, x], [Y_BOT_BAR, Y_BOT_BAR - 0.012], color="#555555",
                    lw=0.6, zorder=3, transform=tr_x)
            ax.text(x, Y_BOT_BAR - 0.018, f"{p / 1e6:g}", ha="center", va="top",
                    fontsize=6.5, color=MB_C, transform=tr_x,
                    rotation_mode="anchor")
            p += step
        ax.text(x0 + w / 2, Y_BOT_BAR - 0.052, name, ha="center", va="top",
                fontsize=9.5, color=qry_color[name], fontweight="bold",
                transform=tr_x, rotation_mode="anchor")

    # ribbons: bezier bands from the bottom edge of the top row to the top
    # edge of the bottom row; minus-strand blocks swap query endpoints
    y_t = Y_TOP_BAR
    y_b = Y_BOT_BAR + BAR_H
    for qn, qs, qe, strand, rn, rs, re_, _na in links:
        ta1, ta2 = ref_frac(rn, rs), ref_frac(rn, re_)
        tb1, tb2 = qry_frac(qn, qs), qry_frac(qn, qe)
        if strand == "-":
            tb1, tb2 = tb2, tb1
        A1, A2 = (ta1, y_t), (ta2, y_t)
        B1, B2 = (tb1, y_b), (tb2, y_b)
        verts = [A1, (A1[0], Y_MID), (B1[0], Y_MID), B1, B2,
                 (B2[0], Y_MID), (A2[0], Y_MID), A2, A1]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                 Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                 Path.CLOSEPOLY]
        ax.add_patch(PathPatch(Path(verts, codes), facecolor=qry_color[qn],
                               edgecolor="none", alpha=args.alpha, zorder=2))

    ax.text(0.5, 0.985, args.title or "Genome synteny", ha="center", va="top",
            fontsize=14, color="#1A1A1A", transform=ax.transAxes,
            rotation_mode="anchor")
    ax.text(0.5, 0.945, f"{len(links)} syntenic blocks | ribbons colored by "
            f"query chromosome | top: reference | bottom: query",
            ha="center", va="top", fontsize=8.5, color="#666666",
            transform=ax.transAxes, rotation_mode="anchor")

    if args.out.lower().endswith(".svg"):
        n_text, _ = svg_export.save_hybrid_svg(fig, args.out)
        print(f"wrote {args.out}  (hybrid SVG: 900-dpi PNG base + {n_text} editable text elements, "
              f"{len(links)} ribbons)")
    else:
        fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight", facecolor="white")
        print(f"wrote {args.out}  ({len(links)} ribbons, "
              f"{len(ref_order)} ref / {len(qry_order)} query sequences)")


if __name__ == "__main__":
    main()
