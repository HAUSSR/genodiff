"""Hybrid SVG export: raster artwork as an embedded high-dpi PNG base, all
text as editable SVG <text> elements.

Alignment guarantee: the figure is rendered at a fixed dpi with text hidden,
and text anchors are captured from the very same matplotlib transform in the
same pixel space (viewBox = PNG pixels). Text therefore sits exactly where
matplotlib placed it; ungrouping in Illustrator/Inkscape keeps positions.
"""
import base64
import io
import xml.sax.saxutils as sx

import matplotlib
import matplotlib.legend as mlegend
import matplotlib.text as mtext

BASE_DPI = 900

HA_MAP = {"left": "start", "center": "middle", "right": "end"}
VA_MAP = {
    "baseline": "auto",
    "bottom": "text-after-edge",
    "center": "central",
    "center_baseline": "central",
    "top": "text-before-edge",
}


def save_hybrid_svg(fig, path: str, dpi: int = BASE_DPI):
    """Render `fig` to a hybrid SVG (raster base + editable text)."""
    fig.set_dpi(dpi)
    fig.canvas.draw()
    W, H = fig.canvas.get_width_height()

    # legend text stays raster (frame + text move together)
    skip = set()
    for leg in fig.findobj(mlegend.Legend):
        skip.update(id(t) for t in leg.findobj(mtext.Text))

    records = []
    for t in fig.findobj(mtext.Text):
        if id(t) in skip or not t.get_visible():
            continue
        s = t.get_text()
        if not s or not s.strip():
            continue
        x, y = t.get_transform().transform(t.get_position())
        fam = t.get_fontfamily()
        fam = fam[0] if isinstance(fam, list) and fam and fam[0] != "sans-serif" else "DejaVu Sans"
        records.append(
            dict(
                s=" ".join(s.split()),  # collapse whitespace/newlines
                x=float(x),
                y=float(y),
                rot=float(t.get_rotation()) % 360,
                ha=t.get_ha(),
                va=t.get_va(),
                fs=float(t.get_size()),
                color=matplotlib.colors.to_hex(t.get_color()),
                fam=fam,
                weight=t.get_weight(),
                style=t.get_fontstyle(),
            )
        )

    # raster pass: hide every text artist, render at the same dpi, no tight
    # bbox (the canvas is the coordinate system and must not change)
    all_text = list(fig.findobj(mtext.Text))
    visible = [t.get_visible() for t in all_text]
    for t in all_text:
        t.set_visible(False)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, facecolor=fig.get_facecolor())
    for t, v in zip(all_text, visible):
        t.set_visible(v)
    png_b64 = base64.b64encode(buf.getvalue()).decode("ascii")

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" '
        f'xmlns:xlink="http://www.w3.org/1999/xlink" '
        f'width="{W / dpi:.4f}in" height="{H / dpi:.4f}in" '
        f'viewBox="0 0 {W} {H}">',
        f'<image x="0" y="0" width="{W}" height="{H}" '
        f'xlink:href="data:image/png;base64,{png_b64}"/>',
    ]
    for r in records:
        fs_px = r["fs"] * dpi / 72.0
        y_svg = H - r["y"]  # SVG y is top-down, matplotlib display y is up
        attrs = [
            f'x="{r["x"]:.2f}"',
            f'y="{y_svg:.2f}"',
            f'font-family="{r["fam"]}"',
            f'font-size="{fs_px:.2f}"',
            f'fill="{r["color"]}"',
            f'text-anchor="{HA_MAP.get(r["ha"], "start")}"',
            f'dominant-baseline="{VA_MAP.get(r["va"], "auto")}"',
        ]
        if r["weight"] and str(r["weight"]).lower() not in ("normal", "0"):
            attrs.append(f'font-weight="{sx.escape(str(r["weight"]))}"')
        if str(r["style"]) == "italic":
            attrs.append('font-style="italic"')
        # positive rotation in matplotlib is CCW (y-up); in SVG y-down it
        # becomes negative. rotation_mode="anchor" in the plots keeps the
        # anchor point identical.
        if r["rot"]:
            attrs.append(f'transform="rotate({-r["rot"]:.2f} {r["x"]:.2f} {y_svg:.2f})"')
        parts.append(f'<text {" ".join(attrs)}>{sx.escape(r["s"])}</text>')
    parts.append("</svg>")

    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(parts))
    return len(records), len(buf.getvalue())
