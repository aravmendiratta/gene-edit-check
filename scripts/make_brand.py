"""Write the brand assets: assets/logo.svg, assets/banner(.-dark).svg, assets/social-card.png.

The mark itself is built in `gene_edit_check/ui/brand.py`. The look and the reasoning are in
docs/DESIGN.md.

    python scripts/make_brand.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch

from gene_edit_check.ui import brand, plots

OUT = Path("assets")
SERIF = "'Iowan Old Style', 'Palatino Linotype', 'Source Serif 4', Georgia, serif"
SANS = "Inter, 'SF Pro Text', 'Segoe UI', 'Helvetica Neue', Arial, sans-serif"

SCHEMES = {
    "": {"panel": "#f4eee0", "ink": "#1b2622", "mute": "#5b655f", "accent": "#1f7a68", "blob": "#1f7a68", "blob2": "#ff7a4d"},
    "-dark": {"panel": "#14201c", "ink": "#efe9dc", "mute": "#a9b0a6", "accent": "#4fbf9f", "blob": "#4fbf9f", "blob2": "#ff7a4d"},
}


def write_logo() -> None:
    (OUT / "logo.svg").write_text(brand.icon_svg(256, uid="logo") + "\n", encoding="utf-8")


def write_banner() -> None:
    for suffix, c in SCHEMES.items():
        svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 240" width="1000" height="240" role="img" aria-label="gene-edit-check: check a gene edit on a computer before the lab">
  <rect width="1000" height="240" rx="28" fill="{c['panel']}"/>
  <g opacity="0.10"><circle cx="905" cy="46" r="86" fill="{c['blob']}"/><circle cx="800" cy="205" r="52" fill="{c['blob2']}"/><circle cx="960" cy="190" r="34" fill="{c['blob']}"/></g>
  <g opacity="0.16" fill="none" stroke="{c['blob']}" stroke-width="1.4"><circle cx="905" cy="46" r="86"/><circle cx="800" cy="205" r="52"/><circle cx="960" cy="190" r="34"/></g>
  {brand.inner_svg(52, 56, 128, "bn" + suffix.strip("-"))}
  <text x="216" y="124" font-family="{SERIF}" font-size="60" font-weight="500" fill="{c['ink']}" letter-spacing="-1.2">gene-edit-<tspan fill="{c['accent']}" font-style="italic">check</tspan></text>
  <text x="218" y="168" font-family="{SANS}" font-size="20" fill="{c['mute']}">Check a gene edit on a computer, before anyone spends weeks on it in a lab.</text>
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
    ax.add_patch(plt.Rectangle((0, 0), 1280, 640, fc="#f4eee0", ec="none"))
    for (cx, cy, r, col, a) in [(1010, 150, 230, "#1f7a68", 0.09), (1160, 520, 130, "#ff7a4d", 0.10), (760, 590, 80, "#1f7a68", 0.08)]:
        ax.add_patch(Circle((cx, cy), r, fc=col, ec="none", alpha=a))
        ax.add_patch(Circle((cx, cy), r, fc="none", ec=col, lw=1.2, alpha=min(1, a * 2.2)))
    brand.draw_mark(ax, 840, 120, 360)
    serif = ["Iowan Old Style", "Palatino Linotype", "Source Serif 4", "Georgia", "DejaVu Serif"]
    ax.text(96, 120, "Computational biology", fontsize=15, color="#5b655f", va="center")
    first = ax.text(96, 300, "gene-edit-", fontsize=62, color="#1b2622", va="center", family=serif)
    fig.canvas.draw()
    end = first.get_window_extent().x1
    ax.text(end - 2, 300, "check", fontsize=62, color="#1f7a68", va="center", family=serif, style="italic")
    ax.text(96, 392, "Check a gene edit on a computer,\nbefore anyone spends weeks on it in a lab.", fontsize=23,
            color="#3d4742", va="top", linespacing=1.45)
    ax.text(96, 556, "Built on AlphaGenome 2  ·  Research use only", fontsize=14, color="#7b837b", va="center")
    fig.savefig(OUT / "social-card.png", dpi=100, facecolor="#f4eee0")
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    for old in OUT.glob("logo-dark.svg"):
        old.unlink()
    write_logo()
    write_banner()
    write_social_card()
    print("wrote", *(str(p) for p in sorted(OUT.iterdir())))
