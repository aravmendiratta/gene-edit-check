"""A fast, 30-second motion-graphics launch video, drawn entirely in code: no screenshots, no stock footage.

    pip install pycairo imageio imageio-ffmpeg numpy
    PYTHONPATH=src python scripts/make_launch_video.py        # writes assets/launch.mp4

Everything on screen is drawn from the project's saved results (`outputs/*.json`): the 8 x 6 check grid, the
48-point efficacy-versus-methylation scatter and the five lab-screen elements are the real numbers. The sound
is synthesised here (kicks, ticks and risers on every cut), so there is nothing to license.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
import wave
from pathlib import Path

import cairo
import numpy as np

W, H, FPS = 1920, 1080, 30
OUT = Path("outputs")

# A louder palette than the figures: this is a video, not a report.
INK = (0.043, 0.043, 0.071)
PAPER = (0.98, 0.97, 0.94)
YELLOW = (1.0, 0.84, 0.04)
PINK = (1.0, 0.24, 0.51)
CYAN = (0.18, 0.89, 0.90)
BLUE = (0.17, 0.29, 1.0)
WHITE = (1, 1, 1)
GREY = (0.45, 0.46, 0.52)

SCENES = [  # name, start, duration, background
    ("hook", 0.0, 3.0, INK),
    ("flip", 3.0, 3.0, YELLOW),
    ("question", 6.0, 3.0, INK),
    ("numbers", 9.0, 4.0, BLUE),
    ("grid", 13.0, 4.0, INK),
    ("scatter", 17.0, 4.0, YELLOW),
    ("lab", 21.0, 4.5, INK),
    ("end", 25.5, 4.5, INK),
]
TOTAL = SCENES[-1][1] + SCENES[-1][2]


# ---- easing -------------------------------------------------------------------------------
def clamp(x: float) -> float:
    return max(0.0, min(1.0, x))


def out_cubic(x: float) -> float:
    x = clamp(x)
    return 1 - (1 - x) ** 3


def in_out(x: float) -> float:
    x = clamp(x)
    return x * x * (3 - 2 * x)


def out_back(x: float, k: float = 1.9) -> float:
    x = clamp(x)
    return 1 + (k + 1) * (x - 1) ** 3 + k * (x - 1) ** 2


def lerp(a: float, b: float, t: float) -> float:
    return a + (b - a) * t


# ---- drawing helpers ----------------------------------------------------------------------
FONT = "Arial Black"


def set_font(c: cairo.Context, size: float, bold: bool = True, family: str | None = None) -> None:
    c.select_font_face(family or FONT, cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD if bold else cairo.FONT_WEIGHT_NORMAL)
    c.set_font_size(size)


def text(c: cairo.Context, s: str, x: float, y: float, size: float, color: tuple, *, align: str = "center",
         scale: float = 1.0, scale_y: float | None = None, alpha: float = 1.0, rot: float = 0.0,
         family: str | None = None, bold: bool = True) -> float:
    """Draw text with its baseline-centre at (x, y). Returns its width. y is the visual centre of the caps."""
    if alpha <= 0.001 or scale <= 0.001:
        return 0.0
    c.save()
    set_font(c, size, bold, family)
    ext = c.text_extents(s)
    c.translate(x, y)
    c.rotate(rot)
    c.scale(scale, scale if scale_y is None else scale_y)
    ox = {"center": -ext.width / 2 - ext.x_bearing, "left": -ext.x_bearing, "right": -ext.width - ext.x_bearing}[align]
    c.move_to(ox, size * 0.36)
    c.set_source_rgba(*color, alpha)
    c.show_text(s)
    c.restore()
    return ext.width * scale


def slam(t: float, t0: float, dur: float = 0.28) -> tuple[float, float]:
    """(scale, alpha) for a word that slams in at t0."""
    if t < t0:
        return 0.0, 0.0
    p = (t - t0) / dur
    return lerp(1.7, 1.0, out_back(p, 1.2)), clamp(p * 4)


def disc(c: cairo.Context, x: float, y: float, r: float, color: tuple, alpha: float = 1.0) -> None:
    if r <= 0 or alpha <= 0:
        return
    c.set_source_rgba(*color, alpha)
    c.arc(x, y, r, 0, 2 * math.pi)
    c.fill()


def line(c: cairo.Context, x0: float, y0: float, x1: float, y1: float, w: float, color: tuple, alpha: float = 1.0) -> None:
    c.set_source_rgba(*color, alpha)
    c.set_line_width(w)
    c.set_line_cap(cairo.LINE_CAP_ROUND)
    c.move_to(x0, y0)
    c.line_to(x1, y1)
    c.stroke()


BASE_COLORS = {"A": CYAN, "T": YELLOW, "C": PINK, "G": WHITE}
PAIR = {"A": "T", "T": "A", "C": "G", "G": "C"}


def helix(c: cairo.Context, cx: float, cy: float, height: float, amp: float, rot: float, *, alpha: float = 1.0,
          n: int = 26, turns: float = 2.4, grow: float = 1.0, seq: str = "ATCGGCTAACGTTGCAGCTA",
          hot: int | None = None, hot_pulse: float = 0.0, dot: float = 11.0) -> None:
    """A rotating double helix of base pairs. `grow` (0..1) draws it from the top down."""
    items = []
    for i in range(n):
        u = i / (n - 1)
        if u > grow:
            break
        y = cy - height / 2 + height * u
        ph = u * turns * 2 * math.pi + rot
        xa, xb = cx + amp * math.sin(ph), cx - amp * math.sin(ph)
        da, db = math.cos(ph), -math.cos(ph)
        b = seq[i % len(seq)]
        items.append((i, y, xa, xb, da, db, b))
    # rungs first (behind), then strands sorted by depth
    for i, y, xa, xb, da, db, b in items:
        col = BASE_COLORS[b]
        is_hot = hot is not None and i == hot
        w = 7 if is_hot else 4
        line(c, xa, y, (xa + xb) / 2, y, w, col, alpha * 0.9)
        line(c, (xa + xb) / 2, y, xb, y, w, BASE_COLORS[PAIR[b]], alpha * 0.9)
        if is_hot and hot_pulse > 0:
            disc(c, (xa + xb) / 2, y, 30 + 120 * hot_pulse, PINK, alpha * 0.5 * (1 - hot_pulse))
    for k in range(len(items) - 1):
        _, y0, xa0, xb0, da0, db0, _ = items[k]
        _, y1, xa1, xb1, da1, db1, _ = items[k + 1]
        line(c, xa0, y0, xa1, y1, 5, WHITE, alpha * (0.22 + 0.3 * max(da0, 0)))
        line(c, xb0, y0, xb1, y1, 5, WHITE, alpha * (0.22 + 0.3 * max(db0, 0)))
    dots = []
    for i, y, xa, xb, da, db, b in items:
        dots += [(da, xa, y, BASE_COLORS[b]), (db, xb, y, BASE_COLORS[PAIR[b]])]
    for d, x, y, col in sorted(dots, key=lambda q: q[0]):
        r = dot * (0.75 + 0.35 * d)
        disc(c, x, y, r, col, alpha * (0.55 + 0.45 * (d + 1) / 2))


def counter(t: float, t0: float, dur: float, to: float) -> float:
    return to * out_cubic((t - t0) / dur)


# ---- data ----------------------------------------------------------------------------------
def load_data() -> dict:
    from gene_edit_check.ui.plots import _CHECKS, _check_row

    rep = json.loads((OUT / "report_crc-pd1-exhaustion.json").read_text(encoding="utf-8"))
    grid = [_check_row(it["report"]["failure_reasons"]) for it in rep["iterations"]]
    pts = [(abs(it["report"]["tier1_selectivity"]["target_delta_sketch_b"]), it["report"]["tier2_hazards"]["cpg_methylation_mean_abs_shift"])
           for it in rep["iterations"]]
    null = rep["efficacy_null"]
    pts += list(zip(null["abs_lfc_values"], null["cpg_methylation_shift_values"]))
    enh = json.loads((OUT / "enhancer_concordance_1000bp.json").read_text(encoding="utf-8"))
    el = enh["profiles"]["target saliency (exhausted TIL)"]["elements"]
    ranks = [(k, v["offset_bp"], v["percentile_vs_null"]) for k, v in el.items()]
    rho = next(x["spearman_rho"] for x in rep["hazard_coupling"]["couplings"] if x["hazard"] == "cpg_methylation")
    return {"grid": grid, "n_checks": len(_CHECKS), "check_names": [c[1] if isinstance(c, (list, tuple)) else str(c) for c in _CHECKS],
            "pts": pts, "ranks": ranks, "rho": rho}


# ---- scenes --------------------------------------------------------------------------------
def scene_hook(c, t, d):
    helix(c, 1500, 540, 960, 200, t * 1.8, alpha=0.8, grow=out_cubic(t / 1.2), dot=17, n=22)
    for s, t0, col, y in (("ONE", 0.15, WHITE, 300), ("LETTER", 0.75, YELLOW, 540), ("OF DNA.", 1.35, PINK, 780)):
        sc, al = slam(t, t0)
        text(c, s, 120, y, 270, col, align="left", scale=sc, alpha=al)


def scene_flip(c, t, d):
    seq = "ATCGCGTAC"
    n = len(seq)
    step = 150
    x0 = W / 2 - (n - 1) * step / 2
    hot = 4
    for i, ch in enumerate(seq):
        t_in = 0.05 + i * 0.05
        sc = out_back((t - t_in) / 0.3, 2.2)
        shown, col = ch, INK
        sy = 1.0
        if i == hot:
            col = PINK
            p = clamp((t - 1.05) / 0.35)
            if p > 0:
                sy = abs(math.cos(p * math.pi))
                shown = "C" if p < 0.5 else "T"
            if p >= 1.0:
                col = INK
        text(c, shown, x0 + i * step, 420, 200, col, scale=max(sc, 0.0), scale_y=max(sc, 0.0) * max(sy, 0.02))
    # flash ring when the letter flips
    pr = clamp((t - 1.2) / 0.6)
    if 0 < pr < 1:
        c.set_source_rgba(*PINK, 1 - pr)
        c.set_line_width(14 * (1 - pr) + 2)
        c.arc(x0 + hot * step, 420, 80 + 520 * out_cubic(pr), 0, 2 * math.pi)
        c.stroke()
    sc, al = slam(t, 0.35)
    text(c, "CHANGE ONE LETTER,", W / 2, 760, 96, INK, scale=sc, alpha=al)
    sc, al = slam(t, 1.65)
    text(c, "TURN A GENE DOWN.", W / 2, 880, 96, PINK, scale=sc, alpha=al)


def scene_question(c, t, d):
    for s, t0, col, y in (("WILL IT WORK", 0.1, WHITE, 190), ("ONLY IN THE CELLS", 0.5, YELLOW, 330), ("YOU WANT?", 0.9, PINK, 470)):
        sc, al = slam(t, t0)
        text(c, s, 120, y, 130, col, align="left", scale=sc, alpha=al)
    rng = np.random.default_rng(3)
    for k in range(34):
        side = k % 2
        bx = 380 if side == 0 else 1380
        px = bx + rng.normal(0, 110)
        py = 760 + rng.normal(0, 70)
        t_in = 0.5 + k * 0.03
        r = 26 * out_back((t - t_in) / 0.3, 2.0)
        col = PINK if side == 0 else CYAN
        hit = clamp((t - 1.9) / 0.15)
        if side == 1 and hit > 0:
            col = YELLOW
            r *= 1 + 0.35 * math.sin((t - 1.9) * 24)
        disc(c, px, py, r, col, 0.95)
    sc, al = slam(t, 0.7)
    text(c, "TUMOUR", 380, 1000, 44, PINK, scale=sc, alpha=al)
    text(c, "HEALTHY", 1380, 1000, 44, CYAN, scale=sc, alpha=al)
    pr = clamp((t - 1.8) / 0.6)
    if 0 < pr < 1:
        c.set_source_rgba(*WHITE, 1 - pr)
        c.set_line_width(8)
        c.arc(880, 760, 40 + 520 * out_cubic(pr), 0, 2 * math.pi)
        c.stroke()


def scene_numbers(c, t, d):
    rng = np.random.default_rng(5)
    for k in range(180):
        x, y = rng.uniform(0, W), rng.uniform(0, H)
        t_in = 0.1 + (k / 180) * 2.8
        r = 7 * out_back((t - t_in) / 0.25, 2.0)
        disc(c, x + 12 * math.sin(t * 2 + k), y + 12 * math.cos(t * 2 + k * 1.3), r, WHITE, 0.22)
    v = counter(t, 0.15, 0.9, 62)
    text(c, f"{v:.0f}", 120, 330, 400, WHITE, align="left", scale=1.0)
    sc, al = slam(t, 0.3)
    text(c, "PATIENTS", 130 + 560, 330, 110, YELLOW, align="left", scale=sc, alpha=al)
    v2 = counter(t, 1.5, 1.1, 12738)
    text(c, f"{v2:,.0f}", 120, 720, 300, WHITE, align="left")
    sc, al = slam(t, 1.7)
    text(c, "EXHAUSTED T CELLS", 130, 900, 90, YELLOW, align="left", scale=sc, alpha=al)


def scene_grid(c, t, d):
    grid = DATA["grid"]
    rows, cols = len(grid), len(grid[0])
    sc, al = slam(t, 0.05)
    text(c, f"{rows} EDITS", 120, 120, 120, WHITE, align="left", scale=sc, alpha=al)
    sc, al = slam(t, 0.35)
    text(c, f"× {cols} SAFETY CHECKS", 120 + 560, 120, 120, YELLOW, align="left", scale=sc, alpha=al)
    gx0, gy0, dx, dy = 960, 250, 140, 82
    hi = clamp((t - 2.1) / 0.3)
    if hi > 0:
        c.set_source_rgba(*YELLOW, 0.16 * hi)
        c.rectangle(gx0 + 3 * dx - dx / 2 + 4, gy0 - 50, dx - 8, dy * rows + 20)
        c.fill()
    fails = 0
    for r, row in enumerate(grid):
        for k, ok in enumerate(row):
            t_in = 0.5 + r * 0.1 + k * 0.025
            p = out_back((t - t_in) / 0.28, 2.4)
            if p <= 0:
                continue
            x, y = gx0 + k * dx, gy0 + r * dy
            if ok:
                c.set_source_rgba(*CYAN, 0.85)
                c.set_line_width(9 * p)
                c.set_line_cap(cairo.LINE_CAP_ROUND)
                c.move_to(x - 20 * p, y + 2 * p)
                c.line_to(x - 5 * p, y + 18 * p)
                c.line_to(x + 24 * p, y - 18 * p)
                c.stroke()
            else:
                col = YELLOW if k == 3 else PINK
                c.set_source_rgba(*col, 1)
                c.set_line_width(10 * p)
                c.set_line_cap(cairo.LINE_CAP_ROUND)
                c.move_to(x - 20 * p, y - 20 * p)
                c.line_to(x + 20 * p, y + 20 * p)
                c.move_to(x + 20 * p, y - 20 * p)
                c.line_to(x - 20 * p, y + 20 * p)
                c.stroke()
    n_meth = sum(1 for row in grid if not row[3])
    sc, al = slam(t, 2.1)
    text(c, f"{n_meth} OF {rows}", 120, 560, 170, YELLOW, align="left", scale=sc, alpha=al)
    text(c, "FAIL ON", 120, 700, 70, WHITE, align="left", scale=sc, alpha=al)
    text(c, "METHYLATION", 120, 780, 70, WHITE, align="left", scale=sc, alpha=al)
    sc, al = slam(t, 3.0, 0.3)
    text(c, "NONE PASSED.", W / 2, 1000, 120, PINK, scale=sc, alpha=al)


def scene_scatter(c, t, d):
    pts = DATA["pts"]
    sc, al = slam(t, 0.05)
    text(c, "STRONGER EDIT", 120, 120, 100, INK, align="left", scale=sc, alpha=al)
    sc, al = slam(t, 0.4)
    text(c, "= MORE DAMAGE", 120, 235, 100, PINK, align="left", scale=sc, alpha=al)
    ax0, ay0, aw, ah = 300, 940, 1400, 600
    c.set_source_rgba(*INK, 0.9)
    line(c, ax0, ay0, ax0 + aw, ay0, 5, INK)
    line(c, ax0, ay0, ax0, ay0 - ah, 5, INK)
    # Spearman is a rank correlation, so the honest picture is the ranks themselves.
    n = len(pts)
    rx = np.argsort(np.argsort([p[0] for p in pts]))
    ry = np.argsort(np.argsort([p[1] for p in pts]))
    sp = [(ax0 + 40 + rx[i] / (n - 1) * (aw - 80), ay0 - 40 - ry[i] / (n - 1) * (ah - 80)) for i in range(n)]
    for i, (px, py) in enumerate(sp):
        t_in = 0.7 + i * 0.03
        r = 15 * out_back((t - t_in) / 0.3, 2.4)
        disc(c, px, py - 40 * (1 - out_cubic((t - t_in) / 0.4)), r, PINK if i < 8 else INK, 0.9)
    # trend line
    p = clamp((t - 2.3) / 0.7)
    if p > 0:
        a = np.polyfit([q[0] for q in sp], [q[1] for q in sp], 1)
        xa, xb = min(q[0] for q in sp), max(q[0] for q in sp)
        line(c, xa, a[0] * xa + a[1], lerp(xa, xb, out_cubic(p)), a[0] * lerp(xa, xb, out_cubic(p)) + a[1], 12, BLUE)
    v = counter(t, 2.4, 0.9, DATA["rho"] if abs(DATA["rho"]) < 1 else 0.57)
    sc, al = slam(t, 2.4, 0.3)
    text(c, f"ρ = +{v:.2f}", 1480, 215, 120, BLUE, scale=sc, alpha=al)
    text(c, "48 EDITS SCORED", 1480, 305, 50, INK, scale=sc, alpha=al)
    text(c, "stronger edit (rank) →", ax0 + aw - 60, ay0 + 50, 34, INK, bold=False, family="Arial", alpha=0.7)
    text(c, "methylation shift (rank) ↑", ax0 - 20, ay0 - ah - 30, 34, INK, bold=False, family="Arial", alpha=0.7, align="left")


def scene_lab(c, t, d):
    sc, al = slam(t, 0.05)
    text(c, "CHECKED AGAINST", 120, 120, 100, WHITE, align="left", scale=sc, alpha=al)
    sc, al = slam(t, 0.4)
    text(c, "A REAL LAB SCREEN.", 120, 240, 100, CYAN, align="left", scale=sc, alpha=al)
    ranks = DATA["ranks"]
    order = sorted(ranks, key=lambda r: r[1])  # genomic order, left to right
    slots = [(n, o, p) for n, o, p in order if o < 0] + [("PDCD1", 0, None)] + [(n, o, p) for n, o, p in order if o > 0]
    x0, x1, base, maxh = 260, 1660, 760, 380
    line(c, x0 - 40, base, x1 + 40, base, 6, WHITE, 0.9)
    gy = base - maxh * 0.95
    a95 = clamp((t - 1.0) / 0.4)
    line(c, x0 - 40, gy, x1 + 40, gy, 3, PINK, 0.8 * a95)
    text(c, "95th percentile", x1 + 40, gy - 26, 32, PINK, align="right", family="Arial", bold=False, alpha=a95)
    step = (x1 - x0) / (len(slots) - 1)
    for i, (name, off, pct) in enumerate(slots):
        x = x0 + i * step
        lab = "PDCD1" if pct is None else (f"{off / 1000:+.2f}".rstrip("0").rstrip(".") + " kb")
        t_in = 0.9 + i * 0.3
        if pct is None:
            disc(c, x, base, 18, YELLOW, clamp((t - 0.7) / 0.3))
            text(c, "PDCD1", x, base + 60, 40, YELLOW, alpha=clamp((t - 0.7) / 0.3))
            continue
        h = maxh * pct * out_cubic((t - t_in) / 0.7)
        found = pct >= 0.95
        col = CYAN if found else GREY
        c.set_source_rgba(*col, 1)
        c.rectangle(x - 46, base - h, 92, h)
        c.fill()
        text(c, lab, x, base + 60, 40, WHITE, alpha=clamp((t - t_in) / 0.3), family="Arial", bold=False)
        if h > 40:
            text(c, "FOUND" if found else "MISSED", x, base - h - 34, 30, CYAN if found else PINK, alpha=clamp((t - t_in - 0.5) / 0.3))
    sc, al = slam(t, 3.1, 0.3)
    text(c, "IT MISSES THE +15 KB ENHANCERS.", W / 2, 960, 84, PINK, scale=sc, alpha=al)


def scene_end(c, t, d):
    helix(c, 470, 540, 940, 165, t * 1.6, grow=out_cubic(t / 1.4), dot=15, n=24)
    sc, al = slam(t, 0.5)
    w = text(c, "gene-edit-", 880, 380, 104, WHITE, align="left", scale=sc, alpha=al)
    text(c, "check", 880 + w, 380, 104, CYAN, align="left", scale=sc, alpha=al)
    sc, al = slam(t, 1.1)
    text(c, "Check a gene edit on a computer", 880, 510, 46, WHITE, align="left", scale=sc, alpha=al, family="Arial")
    text(c, "before anyone spends weeks in a lab.", 880, 575, 46, WHITE, align="left", scale=sc, alpha=al, family="Arial")
    s = "pip install gene-edit-check"
    n = int(len(s) * clamp((t - 1.9) / 0.8))
    text(c, s[:n] + ("_" if int(t * 3) % 2 == 0 else ""), 880, 740, 54, YELLOW, align="left", family="Consolas")
    sc, al = slam(t, 2.9)
    text(c, "github.com/aravmendiratta/gene-edit-check", 880, 840, 38, GREY, align="left", scale=sc, alpha=al, family="Arial", bold=False)


SCENE_FN = {"hook": scene_hook, "flip": scene_flip, "question": scene_question, "numbers": scene_numbers,
            "grid": scene_grid, "scatter": scene_scatter, "lab": scene_lab, "end": scene_end}
DATA: dict = {}


def draw_scene(c: cairo.Context, idx: int, t: float) -> None:
    name, _, dur, bg = SCENES[idx]
    c.set_source_rgb(*bg)
    c.paint()
    c.save()
    z = 1.0 + 0.035 * clamp(t / dur)
    c.translate(W / 2, H / 2)
    c.scale(z, z)
    c.translate(-W / 2, -H / 2)
    SCENE_FN[name](c, t, dur)
    c.restore()


def render(t: float, surf: cairo.ImageSurface) -> None:
    c = cairo.Context(surf)
    idx = max(i for i, s in enumerate(SCENES) if s[1] <= t + 1e-9)
    local = t - SCENES[idx][1]
    if idx > 0 and local < 0.32:
        # wipe: previous scene stays underneath while the new one slides in from the left
        draw_scene(c, idx - 1, SCENES[idx - 1][2] + local)
        c.save()
        c.rectangle(0, 0, W * in_out(local / 0.32), H)
        c.clip()
        draw_scene(c, idx, local)
        c.restore()
    else:
        draw_scene(c, idx, local)
    # soft vignette-free finish: a thin edge flash on each cut
    if idx > 0 and local < 0.12:
        c.set_source_rgba(1, 1, 1, 0.22 * (1 - local / 0.12))
        c.paint()


# ---- sound ---------------------------------------------------------------------------------
def make_audio(path: Path) -> None:
    sr = 44100
    n = int(TOTAL * sr)
    y = np.zeros(n)
    rng = np.random.default_rng(1)

    def add(t0: float, sig: np.ndarray, gain: float = 1.0) -> None:
        i = int(t0 * sr)
        if i >= n:
            return
        m = min(len(sig), n - i)
        y[i : i + m] += sig[:m] * gain

    def kick(dur: float = 0.28) -> np.ndarray:
        tt = np.arange(int(dur * sr)) / sr
        f = 48 + 110 * np.exp(-tt * 28)
        return np.sin(2 * np.pi * np.cumsum(f) / sr) * np.exp(-tt * 9)

    def tick(dur: float = 0.05) -> np.ndarray:
        tt = np.arange(int(dur * sr)) / sr
        return rng.normal(0, 1, len(tt)) * np.exp(-tt * 70) * 0.35

    def riser(dur: float) -> np.ndarray:
        tt = np.arange(int(dur * sr)) / sr
        f = 300 + 3200 * (tt / dur) ** 2
        return np.sin(2 * np.pi * np.cumsum(f) / sr) * (tt / dur) ** 2 * 0.25 + rng.normal(0, 1, len(tt)) * (tt / dur) ** 3 * 0.08

    def blip(freq: float, dur: float = 0.12) -> np.ndarray:
        tt = np.arange(int(dur * sr)) / sr
        return np.sin(2 * np.pi * freq * tt) * np.exp(-tt * 22) * 0.4

    bpm = 120
    beat = 60 / bpm
    k = kick()
    for i in range(int(TOTAL / beat)):
        t0 = i * beat
        add(t0, k, 0.55 if i % 4 == 0 else 0.38)
        add(t0 + beat / 2, tick(), 0.5)
    for _, start, dur, _ in SCENES[1:]:
        add(start - 0.7, riser(0.7), 0.7)
        add(start, kick(0.5), 1.0)
    # a note per slammed word keeps the type feeling struck
    notes = [392, 440, 523, 587, 659]
    j = 0
    for name, start, dur, _ in SCENES:
        for off in (0.15, 0.55, 0.95, 1.45, 2.05):
            if off < dur - 0.3:
                add(start + off, blip(notes[j % len(notes)]), 0.55)
                j += 1
    # simple bass under it all
    root = [55, 55, 65.4, 49]
    for bar in range(int(TOTAL / (beat * 4))):
        f = root[bar % 4]
        tt = np.arange(int(beat * 4 * sr)) / sr
        add(bar * beat * 4, np.sin(2 * np.pi * f * tt) * 0.16 * (0.6 + 0.4 * np.sin(2 * np.pi * tt / beat / 2)), 1.0)
    y = np.tanh(y * 1.1)
    y /= max(1e-6, np.abs(y).max()) / 0.9
    pcm = (y * 32767).astype(np.int16)
    with wave.open(str(path), "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(sr)
        w.writeframes(pcm.tobytes())


def main(out: str = "assets/launch.mp4") -> None:
    import imageio.v2 as imageio
    import imageio_ffmpeg

    DATA.update(load_data())
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    silent = Path(out).with_suffix(".silent.mp4")
    wav = Path(out).with_suffix(".wav")
    surf = cairo.ImageSurface(cairo.FORMAT_ARGB32, W, H)
    writer = imageio.get_writer(str(silent), fps=FPS, codec="libx264", quality=8, macro_block_size=1,
                                ffmpeg_params=["-movflags", "+faststart"])
    total = int(TOTAL * FPS)
    for f in range(total):
        render(f / FPS, surf)
        buf = np.ndarray((H, W, 4), np.uint8, surf.get_data())
        writer.append_data(np.ascontiguousarray(buf[:, :, [2, 1, 0]]))
    writer.close()
    make_audio(wav)
    subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(), "-y", "-loglevel", "error", "-i", str(silent), "-i", str(wav),
                    "-c:v", "copy", "-c:a", "aac", "-b:a", "192k", "-shortest", out], check=True)
    silent.unlink()
    wav.unlink()
    print(f"wrote {out}: {total} frames, {total / FPS:.1f}s, with sound")


def preview(times: list[float], path: str) -> None:
    DATA.update(load_data())
    surf = cairo.ImageSurface(cairo.FORMAT_ARGB32, W, H)
    from PIL import Image

    tiles = []
    for t in times:
        render(t, surf)
        buf = np.ndarray((H, W, 4), np.uint8, surf.get_data())
        tiles.append(Image.fromarray(np.ascontiguousarray(buf[:, :, [2, 1, 0]])).resize((640, 360)))
    cols = 2
    sheet = Image.new("RGB", (640 * cols, 360 * math.ceil(len(tiles) / cols)))
    for i, im in enumerate(tiles):
        sheet.paste(im, ((i % cols) * 640, (i // cols) * 360))
    sheet.save(path)


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--preview":
        preview([float(x) for x in sys.argv[3:]], sys.argv[2])
    else:
        main(*sys.argv[1:2])
