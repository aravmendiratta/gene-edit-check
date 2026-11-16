"""The gene-edit-check mark: a twisting double helix with one coral rung, the edited base.

The geometry is built once here and drawn two ways, as SVG for the README and banner
(`scripts/make_brand.py`) and as matplotlib patches for the figures (`draw_mark`). Everything is in a
100 x 100 box.

The look is meant to feel drawn rather than engineered: ribbons that swell and thin as the strand turns
toward and away from the viewer, a slight lean, and colour that moves along the length instead of a flat fill.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

# Teal and cream with one warm accent. The coral rung is the only warm mark in the picture.
PAPER = "#fbf8f1"
FOREST = "#12453b"
FOREST_DEEP = "#0b2e27"
CREAM = "#f6efe0"
MINT = "#9bd6c0"
TEAL = "#1f7a68"
CORAL = "#ff7a4d"

_TURNS = 1.5
_N = 180
_AMP = 19.0
_LEAN = math.radians(-20)


@dataclass(frozen=True)
class Seg:
    pts: tuple[tuple[float, float], ...]
    depth: float  # -1 (far) .. +1 (near)
    ty: float     # 0 (top) .. 1 (bottom), used for colour drift
    strand: int


@dataclass(frozen=True)
class Rung:
    a: tuple[float, float]
    b: tuple[float, float]
    edited: bool


def _rot(x: float, y: float) -> tuple[float, float]:
    c, s = math.cos(_LEAN), math.sin(_LEAN)
    dx, dy = x - 50, y - 50
    return 50 + dx * c - dy * s, 50 + dx * s + dy * c


def geometry() -> tuple[list[Seg], list[Rung]]:
    total = _TURNS * 2 * math.pi
    segs: list[Seg] = []
    for strand, phase in enumerate((0.0, math.pi)):
        prev = None
        for i in range(_N + 1):
            u = i / _N
            t = u * total
            x = 50 + _AMP * math.sin(t + phase)
            y = 7 + 86 * u
            depth = math.cos(t + phase)
            w = 6.2 + 4.2 * depth  # ribbon swells as it comes toward the viewer
            cur = (x, y, w, depth, u)
            if prev is not None:
                x0, y0, w0, d0, u0 = prev
                quad = (_rot(x0 - w0 / 2, y0), _rot(x - w / 2, y), _rot(x + w / 2, y), _rot(x0 + w0 / 2, y0))
                segs.append(Seg(quad, (d0 + depth) / 2, (u0 + u) / 2, strand))
            prev = cur
    rungs: list[Rung] = []
    ts = [(k + 0.5) * math.pi / 3.2 for k in range(int(total / (math.pi / 3.2)))]
    mid = len(ts) // 2
    for k, t in enumerate(ts):
        u = t / total
        y = 7 + 86 * u
        xa, xb = 50 + _AMP * math.sin(t), 50 + _AMP * math.sin(t + math.pi)
        if abs(xa - xb) < 9:  # skip rungs where the strands cross
            continue
        rungs.append(Rung(_rot(xa, y), _rot(xb, y), edited=(k == mid or (k == mid + 1 and abs(xa - xb) < 9))))
    if not any(r.edited for r in rungs) and rungs:
        j = len(rungs) // 2
        rungs[j] = Rung(rungs[j].a, rungs[j].b, True)
    return segs, rungs


def _mix(c1: str, c2: str, t: float) -> str:
    a = [int(c1[i : i + 2], 16) for i in (1, 3, 5)]
    b = [int(c2[i : i + 2], 16) for i in (1, 3, 5)]
    return "#" + "".join(f"{round(x + (y - x) * t):02x}" for x, y in zip(a, b))


def _strand_color(seg: Seg, front: str, back: str, drift_to: str) -> str:
    base = _mix(front, drift_to, seg.ty)
    near = (seg.depth + 1) / 2  # 0..1
    return _mix(back, base, 0.12 + 0.88 * near)


def _ordered(segs: list[Seg]) -> tuple[list[Seg], list[Seg]]:
    far = sorted((s for s in segs if s.depth < 0), key=lambda s: s.depth)
    near = sorted((s for s in segs if s.depth >= 0), key=lambda s: s.depth)
    return far, near


def icon_svg(size: int = 96, rounded: bool = True, *, bg_top: str = FOREST, bg_bottom: str = FOREST_DEEP,
             uid: str = "g") -> str:
    """Standalone SVG for the app-icon style mark: dark forest tile, cream helix, one coral rung."""
    segs, rungs = geometry()
    far, near = _ordered(segs)
    front, back, drift = CREAM, "#2f7d69", "#a9d9c6"
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100" width="{size}" height="{size}" '
           f'role="img" aria-label="gene-edit-check">',
           f'<defs><linearGradient id="{uid}bg" x1="0" y1="0" x2="0.35" y2="1">'
           f'<stop offset="0" stop-color="{bg_top}"/><stop offset="1" stop-color="{bg_bottom}"/></linearGradient>'
           f'<radialGradient id="{uid}glow" cx="0.5" cy="0.5" r="0.5"><stop offset="0" stop-color="{CORAL}" stop-opacity="0.38"/>'
           f'<stop offset="1" stop-color="{CORAL}" stop-opacity="0"/></radialGradient></defs>']
    out.append(f'<rect width="100" height="100" rx="{22 if rounded else 0}" fill="url(#{uid}bg)"/>')

    def poly(s: Seg) -> str:
        pts = " ".join(f"{x:.2f},{y:.2f}" for x, y in s.pts)
        return f'<polygon points="{pts}" fill="{_strand_color(s, front, back, drift)}" stroke="{_strand_color(s, front, back, drift)}" stroke-width="0.25"/>'

    out += [poly(s) for s in far]
    for r in rungs:
        if r.edited:
            cx, cy = (r.a[0] + r.b[0]) / 2, (r.a[1] + r.b[1]) / 2
            out.append(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="15" fill="url(#{uid}glow)"/>')
            out.append(f'<line x1="{r.a[0]:.2f}" y1="{r.a[1]:.2f}" x2="{r.b[0]:.2f}" y2="{r.b[1]:.2f}" '
                       f'stroke="{CORAL}" stroke-width="3.4" stroke-linecap="round"/>')
            for p in (r.a, r.b):
                out.append(f'<circle cx="{p[0]:.2f}" cy="{p[1]:.2f}" r="3.1" fill="{CORAL}"/>')
        else:
            out.append(f'<line x1="{r.a[0]:.2f}" y1="{r.a[1]:.2f}" x2="{r.b[0]:.2f}" y2="{r.b[1]:.2f}" '
                       f'stroke="{MINT}" stroke-opacity="0.8" stroke-width="2.3" stroke-linecap="round"/>')
    out += [poly(s) for s in near]
    out.append("</svg>")
    return "".join(out)


def inner_svg(x: float, y: float, size: float, uid: str) -> str:
    """The icon as a nested <svg> so it can sit inside a larger SVG (the banner)."""
    s = icon_svg(int(size), uid=uid)
    return s.replace('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 100"',
                     f'<svg x="{x}" y="{y}" viewBox="0 0 100 100"', 1)


def draw_mark(ax: Any, x: float, y: float, size: float, *, tile: bool = True) -> None:
    """Draw the icon into a matplotlib axes whose data units are points or pixels (y grows downward)."""
    from matplotlib.patches import FancyBboxPatch, Polygon

    segs, rungs = geometry()
    far, near = _ordered(segs)
    k = size / 100.0
    front, back, drift = CREAM, "#2f7d69", "#a9d9c6"

    def tr(p: tuple[float, float]) -> tuple[float, float]:
        return x + p[0] * k, y + p[1] * k

    if tile:
        ax.add_patch(FancyBboxPatch((x, y), size, size, boxstyle=f"round,pad=0,rounding_size={size * 0.22}",
                                    fc=FOREST, ec="none", zorder=5))
    for s in far:
        ax.add_patch(Polygon([tr(p) for p in s.pts], closed=True, fc=_strand_color(s, front, back, drift),
                             ec=_strand_color(s, front, back, drift), lw=0.5, zorder=6))
    for r in rungs:
        a, b = tr(r.a), tr(r.b)
        if r.edited:
            ax.plot([a[0], b[0]], [a[1], b[1]], color=CORAL, lw=size * 0.034, solid_capstyle="round", zorder=7)
            ax.scatter([a[0], b[0]], [a[1], b[1]], s=(size * 0.062) ** 2, color=CORAL, zorder=8, linewidths=0)
        else:
            ax.plot([a[0], b[0]], [a[1], b[1]], color=MINT, alpha=0.75, lw=size * 0.02, solid_capstyle="round", zorder=7)
    for s in near:
        ax.add_patch(Polygon([tr(p) for p in s.pts], closed=True, fc=_strand_color(s, front, back, drift),
                             ec=_strand_color(s, front, back, drift), lw=0.5, zorder=9))
