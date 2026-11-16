"""Write the brand assets: assets/logo.svg, assets/banner.svg and assets/social-card.png.

Each SVG comes in a light and a `-dark` variant; the README picks one with a <picture> element. The language behind them is described in docs/DESIGN.md.

    python scripts/make_brand.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from gene_edit_check.ui import plots

OUT = Path("assets")
FONT = "Inter, 'SF Pro Display', 'Segoe UI', 'Helvetica Neue', Arial, sans-serif"

SCHEMES = {
    "": {"ink": "#111111", "mute": "#55534f", "blue": "#2a78d6", "on": "#ffffff"},
    "-dark": {"ink": "#f2f2f0", "mute": "#b4b3ae", "blue": "#5aa0f2", "on": "#0d1117"},
}

# The mark: a rounded square holding a tick. The tick's short stroke is the "base" being edited,
# the long stroke is the verdict.
MARK = """<rect x="{x}" y="{y}" width="{s}" height="{s}" rx="{r}" fill="{blue}"/>
  <path d="M{a} {b} L{c} {d} L{e} {f}" fill="none" stroke="{on}" stroke-width="{w}" stroke-linecap="round" stroke-linejoin="round"/>"""


def mark(x: float, y: float, s: float, c: dict[str, str]) -> str:
    k = s / 100
    return MARK.format(blue=c["blue"], on=c["on"], x=x, y=y, s=s, r=round(s * 0.24, 2), w=round(9 * k, 2),
                       a=round(x + 27 * k, 2), b=round(y + 53 * k, 2),
                       c=round(x + 44 * k, 2), d=round(y + 69 * k, 2),
                       e=round(x + 74 * k, 2), f=round(y + 33 * k, 2))


def write_logo() -> None:
    for suffix, c in SCHEMES.items():
        svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100" width="96" height="96" role="img" '
               f'aria-label="gene-edit-check">\n  {mark(0, 0, 100, c)}\n</svg>\n')
        (OUT / f"logo{suffix}.svg").write_text(svg, encoding="utf-8")


def write_banner() -> None:
    for suffix, c in SCHEMES.items():
        svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 200" width="1000" height="200" role="img" aria-label="gene-edit-check: check a gene edit on a computer before the lab">
  <style>text {{ font-family:{FONT}; }}</style>
  {mark(40, 52, 96, c)}
  <text x="164" y="108" font-size="52" font-weight="600" fill="{c['ink']}" letter-spacing="-1">gene-edit-<tspan fill="{c['blue']}">check</tspan></text>
  <text x="166" y="152" font-size="22" fill="{c['mute']}">Check a gene edit on a computer, before anyone spends weeks on it in a lab.</text>
</svg>
"""
        (OUT / f"banner{suffix}.svg").write_text(svg, encoding="utf-8")


def write_social_card() -> None:
    plots.use_mode("light")
    plots._theme()
    fig = plt.figure(figsize=(12.8, 6.4), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1280)
    ax.set_ylim(640, 0)
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((96, 96), 84, 84, boxstyle="round,pad=0,rounding_size=20", fc=plots.BLUE, ec="none"))
    ax.plot([96 + 23, 96 + 37, 96 + 62], [96 + 45, 96 + 58, 96 + 28], color="white", lw=7, solid_capstyle="round", solid_joinstyle="round")
    first = ax.text(96, 300, "gene-edit-", fontsize=58, fontweight="semibold", color=plots.INK, va="center")
    fig.canvas.draw()
    end = first.get_window_extent().x1  # pixels; the canvas is 100 dpi, 1 data unit = 1 pixel
    ax.text(end - 2, 300, "check", fontsize=58, fontweight="semibold", color=plots.BLUE, va="center")
    ax.text(96, 380, "Check a gene edit on a computer,\nbefore anyone spends weeks on it in a lab.", fontsize=24,
            color=plots.INK_2, va="top", linespacing=1.4)
    ax.text(96, 560, "Computational biology  ·  AlphaGenome 2  ·  Research use only", fontsize=15, color=plots.INK_3, va="center")
    # a quiet motif on the right: the scorecard as a grid of ticks and crosses
    pattern = ["111101", "111101", "110111", "111011", "011101"]
    x0, y0, dx, dy = 800, 150, 62, 62
    for r, row in enumerate(pattern):
        for c, v in enumerate(row):
            cx, cy = x0 + c * dx, y0 + r * dy
            if v == "1":
                ax.plot([cx - 9, cx - 2, cx + 11], [cy + 1, cy + 8, cy - 8], color=plots.NEUTRAL, lw=2.6, solid_capstyle="round", solid_joinstyle="round")
            else:
                ax.plot([cx - 9, cx + 9], [cy - 9, cy + 9], color=plots.ORANGE, lw=2.8, solid_capstyle="round")
                ax.plot([cx - 9, cx + 9], [cy + 9, cy - 9], color=plots.ORANGE, lw=2.8, solid_capstyle="round")
    fig.savefig(OUT / "social-card.png", dpi=100, facecolor=plots.SURFACE)
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    write_logo()
    write_banner()
    write_social_card()
    print("wrote", *(str(p) for p in sorted(OUT.iterdir())))
