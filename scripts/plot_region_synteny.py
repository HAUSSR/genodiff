#!/usr/bin/env python3
"""Region synteny: fragment-level mapping of one reference region into
each query genome (rearrangements, inversions, deletions, translocations).

For every input PAF (a query genome mapped against a common reference),
the chains overlapping the reference region become FRAGMENTS: each fragment
is drawn as its own bar in that genome's row (in genome order, with its
coordinates and strand), connected by one ribbon to the piece of the
reference region it covers. Reading the figure:

  - ribbon crosses (X)          -> fragment inverted relative to the anchor
  - fragment order shuffled     -> rearrangement between the genomes
  - gaps in the anchor coverage -> deletion / diverged candidates
                                   (drawn as colored coverage strips on the
                                   anchor bar's facing edge)
  - fragments on other contigs  -> translocation / paralog candidates

A text report of exactly these findings is printed before the figure is
written.

Rows follow --labels order (N PAFs -> N+1 rows); PAF k connects rows k and
k+1. --orient gives, per PAF, which of its two rows is the PAF's query and
which the reference: "qr" = row k query, row k+1 reference; "rq" flips.
Auto-detected when the two genomes' contig names differ. The region
coordinates are interpreted in the genome named by --region-in.

Usage (Ms region into Mca and Mf; Ms is both PAFs' reference):
    python plot_region_synteny.py Chr3:20725000-20875000 Mca_Ms.paf Mf_Ms.paf \
        --labels Mca,Ms,Mf --orient qr,rq --region-in Ms -o region.svg

Region figures benefit from finer chains than the genome-wide defaults:
    genodiff map ref.gdx query.fa --min-anchors 5 --min-span 200 > q.paf

.svg output is the hybrid export (900-dpi PNG base + editable <text>).
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
from plot_syntenic_linear import PALETTE, natural_key, read_paf

NAME_C, COORD_C = "#1A1A1A", "#666666"
REF_FC, REF_EC = "#D2D2D2", "#7A7A7A"
X0, X1 = 0.125, 0.90          # track bar span (genome names live to the left)
SEG_GAP = 0.012               # gap between fragments of one track


def parse_region(s):
    m = re.fullmatch(r"([\w.\-]+):(\d+)-(\d+)", s.strip())
    if not m:
        raise SystemExit(f"bad region '{s}' (expected chr:start-end, 1-based inclusive)")
    return m.group(1), int(m.group(2)), int(m.group(3))


def fmt_mb(p: float) -> str:
    return f"{p / 1e6:.2f}"


def fmt_span(s0: float, s1: float) -> str:
    """interval label that never collapses to start==end (short fragments
    would read as zero-length at 2 decimals)"""
    a, b = fmt_mb(s0), fmt_mb(s1)
    if a == b and s1 > s0:
        a, b = f"{s0 / 1e6:.3f}", f"{s1 / 1e6:.3f}"
        if a == b:
            a, b = f"{s0 / 1e3:.0f}kb", f"{s1 / 1e3:.0f}kb"
    return f"{a}-{b}"


def nice_step(length: float, max_ticks: int = 6):
    for s in (1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6):
        if length / s <= max_ticks:
            return s
    return 1e7


def merge_intervals(ivs, max_gap):
    out = []
    for s, e in sorted(ivs):
        if out and s - out[-1][1] <= max_gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def split_link(L, side):
    """(contig, start, end) of link L on the given side ('q' or 'r')"""
    return (L[0], L[1], L[2]) if side == "q" else (L[4], L[5], L[6])


def project(links, known_segs, side, merge_gap, min_overlap, min_frag, min_cov_frac):
    """Chains overlapping `known_segs` on side `side` -> (fragments on the
    other genome, kept chains tagged with both segments).

    A chain must overlap the known side by >= min_overlap bp and span >=
    min_frag bp on the other side; fragment clusters are scored by how much
    of the known region they explain, and clusters below min_cov_frac of
    the best cluster's score are dropped as repeat noise."""
    cand = []
    for L in links:
        kc, ks, ke = split_link(L, side)
        oc, os_, oe = split_link(L, "r" if side == "q" else "q")
        if oe - os_ < min_frag:
            continue
        for g in known_segs:
            if g["contig"] != kc:
                continue
            ov = min(ke, g["s1"]) - max(ks, g["s0"])
            if ov >= min_overlap:
                cand.append((L, g, oc, os_, oe, ov))
                break
    per_ctg = {}
    for L, g, oc, os_, oe, ov in cand:
        per_ctg.setdefault(oc, []).append((os_, oe, L[3]))
    frags = []
    for ctg, items in per_ctg.items():
        # a fragment is a COLLINEAR block: same-strand intervals merge only
        # across small gaps AND only while no opposite-strand interval sits
        # between them (an inversion breaks the block)
        items.sort()
        cur = None
        for s, e, st in items:
            if cur is not None and st == cur["strand"] \
                    and s - cur["e"] <= merge_gap:
                cur["e"] = max(cur["e"], e)
            else:
                cur = dict(contig=ctg, s=s, e=e, strand=st, cov=0)
                frags.append(cur)
    for f in frags:
        pad = max(int((f["e"] - f["s"]) * 0.02), 100)
        f["s0"], f["s1"] = max(0, f["s"] - pad), f["e"] + pad
    kept = []
    for L, g, oc, os_, oe, ov in cand:
        mid = (os_ + oe) // 2
        seg = next((f for f in frags
                    if f["contig"] == oc and f["strand"] == L[3]
                    and f["s0"] <= mid <= f["s1"]), None)
        if seg is not None:
            seg["cov"] += ov
            kept.append(dict(L=L, kseg=g, oseg=seg))
    best = max((f["cov"] for f in frags), default=0)
    if best > 0:
        keep = {id(f) for f in frags if f["cov"] >= min_cov_frac * best}
        frags = [f for f in frags if id(f) in keep]
        kept = [c for c in kept if id(c["oseg"]) in keep]
    frags.sort(key=lambda f: (natural_key(f["contig"]), f["s0"]))
    return frags, kept


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("region", help="chr:start-end in --region-in coordinates (1-based)")
    ap.add_argument("pafs", nargs="+", help="PAF per adjacent genome pair, track order")
    ap.add_argument("-o", "--out", default="region_synteny.png")
    ap.add_argument("--labels", required=True,
                    help="comma-separated genome names, N PAFs -> N+1 names")
    ap.add_argument("--orient", default="",
                    help="per-PAF query/ref side of its two rows: qr or rq, "
                         "comma-separated (default: qr for every PAF)")
    ap.add_argument("--region-in", default="",
                    help="genome whose coordinates the region is given in "
                         "(default: first label)")
    ap.add_argument("--min-anchors", type=int, default=3)
    ap.add_argument("--min-overlap", type=int, default=500,
                    help="chain must overlap the region by at least this (bp)")
    ap.add_argument("--min-frag", type=int, default=1000,
                    help="projected fragment must span at least this (bp)")
    ap.add_argument("--min-cov-frac", type=float, default=0.03,
                    help="drop fragment clusters explaining less than this "
                         "fraction of the best cluster (repeat noise); lower "
                         "it to keep small real translocated fragments")
    ap.add_argument("--merge-gap", type=int, default=10_000,
                    help="chains closer than this on the query genome merge "
                         "into one fragment (bp)")
    ap.add_argument("--flank", type=int, default=0,
                    help="extra context on both sides of the region (bp)")
    ap.add_argument("--alpha", type=float, default=0.55)
    ap.add_argument("--title", default="")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    labels = args.labels.split(",")
    if len(labels) != len(args.pafs) + 1:
        raise SystemExit(f"--labels needs {len(args.pafs) + 1} names for "
                         f"{len(args.pafs)} PAFs, got {len(labels)}")
    orients = args.orient.split(",") if args.orient else []
    orients += ["qr"] * (len(args.pafs) - len(orients))
    for o in orients:
        if o not in ("qr", "rq"):
            raise SystemExit(f"bad --orient entry '{o}' (use qr or rq)")
    anchor = args.region_in or labels[0]
    if anchor not in labels:
        raise SystemExit(f"--region-in '{anchor}' is not among --labels")
    ai = labels.index(anchor)

    chrom, rstart, rend = parse_region(args.region)
    r0 = max(0, rstart - 1 - args.flank)   # 0-based half-open
    r1 = rend + args.flank

    pafs = []
    for k, path in enumerate(args.pafs):
        links, ref_len, qry_len = read_paf(path, args.min_anchors)
        pafs.append(dict(path=path, links=links, ref_names=set(ref_len),
                         qry_names=set(qry_len), ribbons=[], sides=None,
                         orient=orients[k]))

    rows = {i: dict(label=labels[i], frags=[]) for i in range(len(labels))}
    rows[ai]["frags"] = [dict(contig=chrom, s0=r0, s1=r1)]

    def fill(k, known_idx):
        unknown_idx = k + 1 if known_idx == k else k
        known_frags = rows[known_idx]["frags"]
        if not known_frags:
            return
        kcs = {g["contig"] for g in known_frags}
        in_q = kcs & pafs[k]["qry_names"]
        in_r = kcs & pafs[k]["ref_names"]
        if in_q and not in_r:
            side = "q"
        elif in_r and not in_q:
            side = "r"
        else:
            side = pafs[k]["orient"][0 if known_idx == k else 1]
        pafs[k]["sides"] = (side, "r" if side == "q" else "q") \
            if known_idx == k else ("r" if side == "q" else "q", side)
        frags, kept = project(pafs[k]["links"], known_frags, side,
                              args.merge_gap, args.min_overlap,
                              args.min_frag, args.min_cov_frac)
        if not frags:
            print(f"  [warn] {pafs[k]['path']}: no chains overlap the region "
                  f"on {labels[known_idx]}", file=sys.stderr)
        rows[unknown_idx]["frags"] = frags
        for c in kept:
            # kseg belongs to row k, k1seg to row k+1, whichever was known
            if known_idx == k:
                pafs[k]["ribbons"].append((c["L"], c["kseg"], c["oseg"]))
            else:
                pafs[k]["ribbons"].append((c["L"], c["oseg"], c["kseg"]))

    for k in range(ai - 1, -1, -1):        # upward from the anchor
        fill(k, k + 1)
    for k in range(ai, len(args.pafs)):    # downward from the anchor
        fill(k, k)

    # ---- text report: per PAF, anchor-side mapping of every fragment ----
    region_bp = r1 - r0
    for k, p in enumerate(pafs):
        a_row = k if ai == k else (k + 1 if ai == k + 1 else None)
        if a_row is None:
            continue  # anchor not part of this pair
        q_row = k + 1 if a_row == k else k
        a_side = p["sides"][0] if a_row == k else p["sides"][1]
        q_side = "r" if a_side == "q" else "q"
        frags = []
        for L, kseg, k1seg in p["ribbons"]:
            ac, as_, ae = split_link(L, a_side)
            qc, qs, qe = split_link(L, q_side)
            as_, ae = max(as_, r0), min(ae, r1)
            frags.append(dict(a=(as_, ae), q=(qc, qs, qe), strand=L[3]))
        frags.sort(key=lambda f: f["a"][0])
        n_blocks = len({id(k1seg if q_row == k + 1 else kseg)
                        for _, kseg, k1seg in p["ribbons"]})
        cov = merge_intervals([(f["a"][0], f["a"][1]) for f in frags], 0)
        cov_bp = sum(e - s for s, e in cov)
        print(f"\n{labels[q_row]} vs {labels[a_row]}  {chrom}:"
              f"{fmt_mb(r0)}-{fmt_mb(r1)} Mb")
        print(f"  blocks: {n_blocks} | chains: {len(frags)} | anchor coverage: "
              f"{cov_bp}/{region_bp} bp ({100 * cov_bp / max(region_bp, 1):.0f}%)")
        for f in frags:
            print(f"    anchor {fmt_span(f['a'][0], f['a'][1])} Mb"
                  f"  <->  {f['q'][0]}:{fmt_span(f['q'][1], f['q'][2])} Mb"
                  f"  {'+' if f['strand'] == '+' else '- (inverted)'}")
        gaps = []
        pos = r0
        for s, e in cov:
            if s - pos >= 1000:
                gaps.append(f"{fmt_span(pos, s)} Mb ({(s - pos) / 1e3:.1f} kb)")
            pos = max(pos, e)
        if r1 - pos >= 1000:
            gaps.append(f"{fmt_span(pos, r1)} Mb ({(r1 - pos) / 1e3:.1f} kb)")
        if gaps:
            print(f"  uncovered in {labels[q_row]} (deletion/diverged candidates): "
                  + "; ".join(gaps))
        inv = [fmt_span(f["a"][0], f["a"][1]) for f in frags
               if f["strand"] == "-"]
        if inv:
            print(f"  inverted fragments (anchor Mb): " + ", ".join(inv))
        main_ctg = max((f["q"][0] for f in frags),
                       key=lambda c: sum(f["a"][1] - f["a"][0]
                                         for f in frags if f["q"][0] == c),
                       default=None)
        other = [f"{f['q'][0]}:{fmt_span(f['q'][1], f['q'][2])}"
                 for f in frags if f["q"][0] != main_ctg]
        if other:
            print(f"  fragments off {main_ctg} (translocation/paralog candidates): "
                  + ", ".join(other))
        order = sorted(range(len(frags)),
                       key=lambda i: (natural_key(frags[i]["q"][0]),
                                      frags[i]["q"][1], frags[i]["q"][2]))
        if frags and order != list(range(len(frags))):
            print(f"  fragment order differs from anchor order "
                  f"(anchor rank -> genome rank): "
                  + ", ".join(str(i + 1) for i in order))

    # ---- figure ----
    n_rows = len(labels)
    fig, ax = plt.subplots(figsize=(14, 1.9 * n_rows + 1.6))
    ax.set_position([0.01, 0.01, 0.985, 0.98])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    tr_x = ax.get_xaxis_transform()

    y_top, y_bot = 0.90, 0.04
    row_h = (y_top - y_bot) / n_rows
    bar_h = min(0.045, 0.55 * row_h)
    row_y = [y_top - (i + 1) * row_h + (row_h - bar_h) / 2 for i in range(n_rows)]
    colors = {i: (REF_EC if i == ai else PALETTE[i % len(PALETTE)])
              for i in range(n_rows)}
    rib_color = {k: colors[k if ai == k + 1 else (k + 1 if ai == k else k + 1)]
                 for k in range(len(pafs))}

    def layout(frags):
        total = sum(f["s1"] - f["s0"] for f in frags) or 1
        avail = X1 - X0 - SEG_GAP * (len(frags) - 1)
        x = X0
        for f in frags:
            w = avail * (f["s1"] - f["s0"]) / total
            f["x0"], f["x1"] = x, x + w
            x += w + SEG_GAP

    def frac(f, pos):
        pos = min(max(pos, f["s0"]), f["s1"])
        return f["x0"] + (pos - f["s0"]) / (f["s1"] - f["s0"]) * (f["x1"] - f["x0"])

    for i in range(n_rows):
        layout(rows[i]["frags"])

    def band(ta1, ta2, tb1, tb2, y_a, y_b, color):
        y_mid = (y_a + y_b) / 2
        A1, A2 = (ta1, y_a), (ta2, y_a)
        B1, B2 = (tb1, y_b), (tb2, y_b)
        verts = [A1, (A1[0], y_mid), (B1[0], y_mid), B1, B2,
                 (B2[0], y_mid), (A2[0], y_mid), A2, A1]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                 Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                 Path.CLOSEPOLY]
        ax.add_patch(PathPatch(Path(verts, codes), facecolor=color,
                               edgecolor="none", alpha=args.alpha, zorder=2))

    # coverage strips on the anchor bar facing each query (gaps = deletions)
    STRIP = 0.011
    for k, p in enumerate(pafs):
        a_row = k if ai == k else (k + 1 if ai == k + 1 else None)
        if a_row is None:
            continue
        a_side = p["sides"][0] if a_row == k else p["sides"][1]
        afrag = rows[ai]["frags"][0]
        y = row_y[ai]
        y_edge = y + bar_h if k < ai else y - STRIP  # facing the query row
        for L, kseg, k1seg in p["ribbons"]:
            _, as_, ae = split_link(L, a_side)
            x1 = frac(afrag, max(as_, r0))
            x2 = frac(afrag, min(ae, r1))
            if x2 > x1:
                ax.add_patch(Rectangle((x1, y_edge), x2 - x1, STRIP,
                                       fc=rib_color[k], ec="none",
                                       transform=tr_x, zorder=4))

    # ribbons: one per chain
    for k, p in enumerate(pafs):
        y_a = row_y[k]
        y_b = row_y[k + 1] + bar_h
        for L, kseg, k1seg in p["ribbons"]:
            ca, a1, a2 = split_link(L, p["sides"][0])
            cb, b1, b2 = split_link(L, p["sides"][1])
            tb1, tb2 = frac(k1seg, b1), frac(k1seg, b2)
            if L[3] == "-":
                tb1, tb2 = tb2, tb1
            band(frac(kseg, a1), frac(kseg, a2), tb1, tb2, y_a, y_b,
                 rib_color[k])

    # fragment strand relative to the anchor (for the +/− label suffix)
    for i in range(n_rows):
        if i == ai:
            continue
        k = i - 1 if ai < i else i
        if not (0 <= k < len(pafs)):
            continue
        for L, kseg, k1seg in pafs[k]["ribbons"]:
            f = k1seg if ai < i else kseg
            f["strand_n"] = f.get("strand_n", 0) + (1 if L[3] == "+" else -1)

    # bars, labels
    for i in range(n_rows):
        y = row_y[i]
        for j, f in enumerate(rows[i]["frags"]):
            ax.add_patch(Rectangle((f["x0"], y), f["x1"] - f["x0"], bar_h,
                                   fc=REF_FC if i == ai else colors[i],
                                   ec=REF_EC if i == ai else "#333333",
                                   lw=0.8, zorder=3, transform=tr_x))
            if i != ai:
                st = f.get("strand_n", 0)
                ax.text((f["x0"] + f["x1"]) / 2,
                        y - 0.030 - 0.024 * (j % 2),
                        f"{f['contig']}:{fmt_span(f['s0'], f['s1'])}"
                        f" ({'+' if st >= 0 else '-'})",
                        ha="center", va="top", fontsize=6, color=COORD_C,
                        transform=tr_x)
            elif len(rows[i]["frags"]) == 1:
                step = nice_step(f["s1"] - f["s0"])
                fmt = (lambda t: f"{t / 1e6:.3f}") if step < 1e4 else \
                      (lambda t: f"{t / 1e6:.2f}")
                t = (f["s0"] // step + 1) * step
                while t < f["s1"] * 0.999:
                    x = frac(f, t)
                    ax.plot([x, x], [y, y - 0.012], color="#555555", lw=0.6,
                            zorder=3, transform=tr_x)
                    ax.text(x, y - 0.016, fmt(t), ha="center",
                            va="top", fontsize=6.5, color=COORD_C,
                            transform=tr_x)
                    t += step
        ax.text(0.012, y + bar_h / 2, labels[i], ha="left", va="center",
                fontsize=10, color=NAME_C if i == ai else colors[i],
                fontweight="bold", transform=tr_x)
        if i == ai:
            ax.text(X1 + 0.006, y + bar_h / 2,
                    f"{chrom}:{fmt_mb(r0)}-{fmt_mb(r1)} Mb (region)",
                    ha="left", va="center", fontsize=7.5, color=COORD_C,
                    transform=tr_x)
        if not rows[i]["frags"]:
            ax.text((X0 + X1) / 2, y + bar_h / 2, "no syntenic fragment",
                    ha="center", va="center", fontsize=9, color="#999999",
                    style="italic", transform=tr_x)

    region_str = f"{chrom}:{fmt_mb(r0)}-{fmt_mb(r1)} Mb in {anchor}"
    ax.text(0.5, 0.985, args.title or f"Region synteny {region_str}",
            ha="center", va="top", fontsize=14, color=NAME_C,
            transform=ax.transAxes)
    ax.text(0.5, 0.955, "one ribbon per syntenic fragment | crossing = inversion | "
            "color strip gaps on the anchor bar = deletion/diverged | shuffled "
            "fragment order = rearrangement | other-contig fragments = "
            "translocation candidates", ha="center", va="top", fontsize=8,
            color="#666666", transform=ax.transAxes)

    if args.out.lower().endswith(".svg"):
        n_text, _ = svg_export.save_hybrid_svg(fig, args.out)
        print(f"\nwrote {args.out}  (hybrid SVG: 900-dpi PNG base + "
              f"{n_text} editable text elements)")
    else:
        fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight", facecolor="white")
        print(f"\nwrote {args.out}  ({sum(len(p['ribbons']) for p in pafs)} "
              f"fragment chains across {n_rows} genomes)")


if __name__ == "__main__":
    main()
