#!/usr/bin/env python3
"""Syntenic pairs dot plot from a PAF file (genodiff map / minimap2 compatible).

Each alignment block is one line segment in the query x reference plane:
plus-strand blocks in blue, minus-strand blocks (inversions) in red. The
reference genome runs along the bottom x axis, the query along the left y
axis (both with chromosome bands, names and Mb scales).

Output .png/.pdf renders a plain figure; output .svg renders a hybrid SVG:
900-dpi PNG base with all text as separately editable, exactly anchored
<text> elements.

Usage:
    python plot_syntenic_dotplot.py input.paf -o dotplot.png
    python plot_syntenic_dotplot.py input.paf -o dotplot.svg --top 8
"""
import argparse
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svg_export

PLUS_C = "#2166AC"    # strong blue
MINUS_C = "#B2182B"   # strong red
BAND_C = "#F0F0F0"
SEP_C = "#A8A8A8"
NAME_C = "#1A1A1A"
MB_C = "#666666"


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
                qlen, rlen = int(f[1]), int(f[6])
                na = int(f[9]) if len(f) > 9 and f[9].lstrip("-").isdigit() else 1
            except ValueError:
                continue
            if na < min_anchors:
                continue
            links.append((f[0], qs, qe, f[4], f[5], rs, re_, na))
            qry_len[f[0]] = qlen
            ref_len[f[5]] = rlen
    return links, ref_len, qry_len


def layout(lengths: dict, order: list, equal: bool):
    off, offsets = 0, {}
    size = max(lengths.values()) if equal else 0
    for name in order:
        offsets[name] = off
        off += size if equal else lengths[name]
    return offsets, off


def nice_mb_step(total_bp: float, max_ticks: int = 10):
    mb = 1e6
    total = total_bp / mb
    for step in (0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500):
        if total / step <= max_ticks:
            return step * mb
    return 1000 * mb


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paf", help="input PAF (genodiff map output)")
    ap.add_argument("-o", "--out", default="syntenic_dotplot.png")
    ap.add_argument("--min-anchors", type=int, default=25,
                    help="ignore blocks with fewer anchors (default 25)")
    ap.add_argument("--top", type=int, default=0,
                    help="keep only the N longest query sequences (0 = all)")
    ap.add_argument("--equal", action="store_true",
                    help="give every chromosome an equal-width band (nice with many scaffolds)")
    ap.add_argument("--alpha", type=float, default=0.75)
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

    roff, ref_total = layout(ref_len, ref_order, args.equal)
    qoff, qry_total = layout(qry_len, qry_order, args.equal)
    size = lambda d, n: max(d[n], 1) if args.equal else d[n]
    mb = 1e6

    fig, ax = plt.subplots(figsize=(11, 11))
    fig.subplots_adjust(left=0.155, right=0.975, top=0.925, bottom=0.155)

    # alternating chromosome bands + separators
    for i, name in enumerate(ref_order):
        x0, w = roff[name], size(ref_len, name)
        if i % 2 == 0:
            ax.add_patch(plt.Rectangle((x0, 0), w, qry_total, fc=BAND_C, ec="none", zorder=0))
        ax.plot([x0, x0], [0, qry_total], color=SEP_C, lw=0.7, zorder=1)
    for i, name in enumerate(qry_order):
        y0, h = qoff[name], size(qry_len, name)
        if i % 2 == 0:
            ax.add_patch(plt.Rectangle((0, y0), ref_total, h, fc=BAND_C, ec="none", zorder=0))
        ax.plot([0, ref_total], [y0, y0], color=SEP_C, lw=0.7, zorder=1)

    # alignment blocks
    segs, cols, lws = [], [], []
    for qn, qs, qe, strand, rn, rs, re_, na in links:
        if qn not in qoff or rn not in roff:
            continue
        x0, x1 = roff[rn] + rs, roff[rn] + re_
        y0, y1 = (qoff[qn] + qs, qoff[qn] + qe) if strand == "+" else (qoff[qn] + qe, qoff[qn] + qs)
        segs.append(((x0, y0), (x1, y1)))
        cols.append(PLUS_C if strand == "+" else MINUS_C)
        lws.append(0.7 + min(2.4, (re_ - rs) ** 0.5 / 80.0))
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=lws,
                                     alpha=args.alpha, zorder=2))

    ax.set_xlim(0, ref_total)
    ax.set_ylim(0, qry_total)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    # Mb scales directly under the plot (x) and left of it (y)
    tr_x = ax.get_xaxis_transform()   # x in data, y in axes fraction
    tr_y = ax.get_yaxis_transform()   # x in axes fraction, y in data
    step_r = nice_mb_step(ref_total)
    x = 0.0
    while x <= ref_total:
        ax.text(x, -0.012, f"{x / mb:g}", transform=tr_x, ha="center", va="top",
                fontsize=8, color=MB_C, rotation_mode="anchor")
        x += step_r
    step_q = nice_mb_step(qry_total)
    y = 0.0
    while y <= qry_total:
        ax.text(-0.012, y, f"{y / mb:g}", transform=tr_y, ha="right", va="center",
                fontsize=8, color=MB_C, rotation_mode="anchor")
        y += step_q

    # chromosome names: bottom row and left column, centred per segment
    rot_x = 45 if len(ref_order) > 12 else 0
    for name in ref_order:
        cx = roff[name] + size(ref_len, name) / 2
        ax.text(cx, -0.045, name, transform=tr_x, ha="center", va="top",
                fontsize=9.5, color=NAME_C, fontweight="bold",
                rotation=rot_x, rotation_mode="anchor")
    for name in qry_order:
        cy = qoff[name] + size(qry_len, name) / 2
        ax.text(-0.045, cy, name, transform=tr_y, ha="right", va="center",
                fontsize=9.5, color=NAME_C, fontweight="bold", rotation_mode="anchor")

    ax.set_xlabel(f"Reference (Mb) | {len(ref_order)} sequences", fontsize=11, labelpad=48)
    ax.set_ylabel(f"Query (Mb) | {len(qry_order)} sequences", fontsize=11, labelpad=48)
    ax.set_title(args.title or "Syntenic dot plot", fontsize=14, pad=16)

    ax.legend(
        handles=[
            Line2D([], [], color=PLUS_C, lw=2.4, label="plus-strand block"),
            Line2D([], [], color=MINUS_C, lw=2.4, label="minus-strand block (inverted)"),
        ],
        loc="upper left", frameon=True, framealpha=0.95, edgecolor="#CCCCCC", fontsize=9,
    )

    if args.out.lower().endswith(".svg"):
        n_text, n_png = svg_export.save_hybrid_svg(fig, args.out)
        print(f"wrote {args.out}  (hybrid SVG: 900-dpi PNG base + {n_text} editable text elements, "
              f"{len(segs)} blocks)")
    else:
        fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight", facecolor="white")
        print(f"wrote {args.out}  ({len(segs)} blocks, "
              f"{len(ref_order)} ref / {len(qry_order)} query sequences)")


if __name__ == "__main__":
    main()
