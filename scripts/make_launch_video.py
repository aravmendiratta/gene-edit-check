"""A 90-second motion-graphics explainer for gene-edit-check, drawn entirely in code.

    pip install pycairo imageio-ffmpeg numpy scipy
    PYTHONPATH=src python scripts/make_launch_video.py        # writes assets/launch.mp4 (1280x720)

Structure: the tool is introduced first, then each step of the pipeline gets its own heading and two to four
animated beats under it, then the result, then how to use it. No screenshots and no stock footage. The numbers
(the random-edit null, the candidates, the patient and cell counts) are read from the saved results in
`outputs/`. The music is synthesised here (112 BPM, A minor), so nothing needs licensing. Picture and sound are
muxed in one ffmpeg pass with the index at the front, so the file streams without stalling.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import cairo
import numpy as np

W, H = 1920, 1080            # the layout grid everything is drawn on
OW, OH, FPS = 1280, 720, 30  # what gets rendered
BPM = 112
BEAT = 60 / BPM
BAR = BEAT * 4
OUT = Path("outputs")

# Deep forest, warm paper, soft teal, one muted coral for anything that fails.
FOREST = (0.055, 0.130, 0.115)
FOREST_2 = (0.085, 0.175, 0.155)
PAPER = (0.957, 0.933, 0.875)
INK = (0.090, 0.150, 0.130)
TEAL = (0.180, 0.620, 0.530)
MINT = (0.610, 0.840, 0.750)
CORAL = (0.929, 0.450, 0.310)
GOLD = (0.910, 0.770, 0.420)
LILAC = (0.650, 0.600, 0.820)
MUTE_D = (0.600, 0.690, 0.650)
MUTE_L = (0.380, 0.450, 0.420)
CARD = (0.985, 0.972, 0.935)
TERM = (0.043, 0.105, 0.092)

SCENES = [  # name, bars, dark background?
    ("intro", 4, True), ("find", 5, False), ("predict", 6, True), ("stress", 6, False),
    ("judge", 7, True), ("stats", 3, False), ("use", 8, True), ("end", 3, True),
]
STARTS = np.cumsum([0.0] + [b * BAR for _, b, _ in SCENES])
TOTAL = float(STARTS[-1])
DATA: dict = {}
SERIF, SANS, MONO = "Georgia", "Segoe UI", "Consolas"


# ---- easing ------------------------------------------------------------------------------------
def clamp(x: float) -> float:
    return max(0.0, min(1.0, x))


def out_cubic(x: float) -> float:
    return 1 - (1 - clamp(x)) ** 3


def in_out(x: float) -> float:
    x = clamp(x)
    return x * x * (3 - 2 * x)


def out_back(x: float, k: float = 1.6) -> float:
    x = clamp(x)
    return 1 + (k + 1) * (x - 1) ** 3 + k * (x - 1) ** 2


def lerp(a: float, b: float, t: float) -> float:
    return a + (b - a) * t


def fade(t: float, t0: float, dur: float = 0.7) -> float:
    return in_out((t - t0) / dur)


def rise(t: float, t0: float, dur: float = 0.8, dist: float = 26.0) -> tuple[float, float]:
    p = out_cubic((t - t0) / dur)
    return p, (1 - p) * dist


# ---- primitives --------------------------------------------------------------------------------
def text(c, s, x, y, size, color, *, align="center", alpha=1.0, family=SANS, bold=True, dy=0.0, rot=0.0, fit=True):
    if alpha <= 0.004:
        return 0.0
    c.save()
    c.select_font_face(family, cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD if bold else cairo.FONT_WEIGHT_NORMAL)
    c.set_font_size(size)
    ext = c.text_extents(s)
    sc = 1.0
    if fit and ext.width > 0:
        avail = {"left": W - x - 90, "right": x - 90}.get(align, 2 * min(x, W - x) - 120)
        if ext.width > avail:
            sc = avail / ext.width
    c.translate(x, y + dy)
    c.rotate(rot)
    c.scale(sc, sc)
    ox = {"center": -ext.width / 2 - ext.x_bearing, "left": -ext.x_bearing, "right": -ext.width - ext.x_bearing}[align]
    c.move_to(ox, size * 0.34)
    c.set_source_rgba(*color, alpha)
    c.show_text(s)
    c.restore()
    return ext.width * sc


def disc(c, x, y, r, color, alpha=1.0):
    if r <= 0 or alpha <= 0:
        return
    c.new_path()
    c.set_source_rgba(*color, alpha)
    c.arc(x, y, r, 0, 2 * math.pi)
    c.fill()


def ring(c, x, y, r, color, w, alpha=1.0):
    if r <= 0 or alpha <= 0:
        return
    c.new_path()
    c.set_source_rgba(*color, alpha)
    c.set_line_width(w)
    c.arc(x, y, r, 0, 2 * math.pi)
    c.stroke()


def line(c, x0, y0, x1, y1, w, color, alpha=1.0, dash=None):
    c.new_path()
    c.set_source_rgba(*color, alpha)
    c.set_line_width(w)
    c.set_line_cap(cairo.LINE_CAP_ROUND)
    if dash:
        c.set_dash(dash)
    c.move_to(x0, y0)
    c.line_to(x1, y1)
    c.stroke()
    c.set_dash([])


def rrect(c, x, y, w, h, r, color, alpha=1.0, stroke=0.0):
    if w <= 0 or h <= 0 or alpha <= 0:
        return
    r = min(r, w / 2, h / 2)
    c.new_path()
    c.arc(x + w - r, y + r, r, -math.pi / 2, 0)
    c.arc(x + w - r, y + h - r, r, 0, math.pi / 2)
    c.arc(x + r, y + h - r, r, math.pi / 2, math.pi)
    c.arc(x + r, y + r, r, math.pi, 3 * math.pi / 2)
    c.close_path()
    c.set_source_rgba(*color, alpha)
    if stroke:
        c.set_line_width(stroke)
        c.stroke()
    else:
        c.fill()


def tick(c, x, y, p, color, w=9):
    if p <= 0:
        return
    c.set_source_rgba(*color, clamp(p * 1.5))
    c.set_line_width(w)
    c.set_line_cap(cairo.LINE_CAP_ROUND)
    c.set_line_join(cairo.LINE_JOIN_ROUND)
    c.move_to(x - 20 * p, y + 2 * p)
    c.line_to(x - 5 * p, y + 17 * p)
    c.line_to(x + 22 * p, y - 17 * p)
    c.stroke()


def cross(c, x, y, p, color, w=9):
    if p <= 0:
        return
    c.set_source_rgba(*color, clamp(p * 1.5))
    c.set_line_width(w)
    c.set_line_cap(cairo.LINE_CAP_ROUND)
    c.move_to(x - 17 * p, y - 17 * p)
    c.line_to(x + 17 * p, y + 17 * p)
    c.move_to(x + 17 * p, y - 17 * p)
    c.line_to(x - 17 * p, y + 17 * p)
    c.stroke()


def base_color(ch: str, dark: bool = True):
    return {"A": MINT if dark else TEAL, "T": GOLD, "C": CORAL, "G": PAPER if dark else INK}[ch]


PAIR = {"A": "T", "T": "A", "C": "G", "G": "C"}


def helix(c, cx, cy, height, amp, rot, *, alpha=1.0, n=26, turns=2.2, grow=1.0,
          seq="ATCGGCTAACGTTGCAGCTA", dot=11.0, dark=True):
    items = []
    for i in range(n):
        u = i / (n - 1)
        if u > grow:
            break
        y = cy - height / 2 + height * u
        ph = u * turns * 2 * math.pi + rot
        items.append((y, cx + amp * math.sin(ph), cx - amp * math.sin(ph), math.cos(ph), seq[i % len(seq)]))
    for y, xa, xb, d, b in items:
        line(c, xa, y, (xa + xb) / 2, y, 4, base_color(b, dark), alpha * 0.8)
        line(c, (xa + xb) / 2, y, xb, y, 4, base_color(PAIR[b], dark), alpha * 0.8)
    fg = PAPER if dark else INK
    for k in range(len(items) - 1):
        y0, xa0, xb0, d0, _ = items[k]
        y1, xa1, xb1, _, _ = items[k + 1]
        line(c, xa0, y0, xa1, y1, 5, fg, alpha * (0.15 + 0.25 * max(d0, 0)))
        line(c, xb0, y0, xb1, y1, 5, fg, alpha * (0.15 + 0.25 * max(-d0, 0)))
    dots = []
    for y, xa, xb, d, b in items:
        dots += [(d, xa, y, base_color(b, dark)), (-d, xb, y, base_color(PAIR[b], dark))]
    for d, x, y, col in sorted(dots, key=lambda q: q[0]):
        disc(c, x, y, dot * (0.75 + 0.3 * d), col, alpha * (0.5 + 0.5 * (d + 1) / 2))


def strand(c, x0, x1, y, seq, dark, *, size=34, alpha=1.0, highlight=None, hl_color=CORAL, reveal=1.0):
    """A row of DNA letters. `highlight` is a set of indices drawn in `hl_color`."""
    n = len(seq)
    step = (x1 - x0) / max(n - 1, 1)
    for i, ch in enumerate(seq):
        if i / n > reveal:
            break
        hl = highlight is not None and i in highlight
        col = hl_color if hl else (base_color(ch, dark) if dark else base_color(ch, False))
        text(c, ch, x0 + i * step, y, size, col, alpha=alpha * (1.0 if hl else 0.85), family=MONO, fit=False)
    return step


def cell(c, x, y, r, col, *, alpha=1.0, fill=1.0, seed=0, t=0.0, n=26, dark=True):
    """A cell outline with `fill` of its expression dots still inside."""
    disc(c, x, y, r, col, alpha * 0.12)
    ring(c, x, y, r, col, 5, alpha)
    rng = np.random.default_rng(seed)
    for j in range(n):
        ang, rr = rng.uniform(0, 2 * math.pi), (r - 24) * math.sqrt(rng.uniform(0.05, 1))
        keep = rng.uniform(0, 1) < fill
        if not keep:
            continue
        disc(c, x + rr * math.cos(ang) + 4 * math.sin(t * 1.1 + j), y + rr * math.sin(ang) + 4 * math.cos(t * 0.9 + j),
             7, PAPER if dark else INK, alpha * 0.8)


def counter(t, t0, dur, to):
    return to * out_cubic((t - t0) / dur)


# ---- scene scaffolding -------------------------------------------------------------------------
def header(c, t, dur, n, title, dark, progress=True):
    """A step heading that stays put for the whole scene, with a thin progress line under it."""
    fg = PAPER if dark else INK
    a, dy = rise(t, 0.15, 0.8)
    if n:
        disc(c, 150, 142 + dy, 32, CORAL, a)
        text(c, str(n), 150, 142 + dy, 34, PAPER, alpha=a)
    x0 = 212 if n else 150
    text(c, title, x0, 142, 72, fg, align="left", alpha=a, family=SERIF, bold=False, dy=dy)
    if progress:
        y = 196
        line(c, x0, y, W - 150, y, 3, fg, a * 0.14)
        p = clamp(t / dur)
        if p > 0:
            line(c, x0, y, x0 + (W - 150 - x0) * p, y, 3, CORAL, a * 0.9)


def stage(c, t, beats, fade_dur=0.65):
    """Crossfade a list of (duration, draw_fn, caption) beats. Returns the caption to show and its alpha."""
    starts = [0.0]
    for d, _, _ in beats:
        starts.append(starts[-1] + d)
    best = ("", 0.0)
    for i, (d, fn, cap) in enumerate(beats):
        t_in = starts[i] - (fade_dur if i else 0.0)
        t_out = starts[i + 1]
        if t < t_in or (i < len(beats) - 1 and t >= t_out):
            continue
        a = 1.0
        if i:
            a = min(a, in_out((t - t_in) / fade_dur))
        if i < len(beats) - 1:
            a = min(a, 1 - in_out((t - (t_out - fade_dur)) / fade_dur))
        if a <= 0.004:
            continue
        c.push_group()
        fn(c, t - t_in)
        c.pop_group_to_source()
        c.paint_with_alpha(a)
        if a > best[1]:
            best = (cap, a)
    return best


def caption(c, cap, alpha, dark, y=990):
    if cap:
        text(c, cap, W / 2, y, 40, MUTE_D if dark else MUTE_L, alpha=alpha, bold=False)


# ---- data --------------------------------------------------------------------------------------
def load_data() -> dict:
    rep = json.loads((OUT / "report_crc-pd1-exhaustion.json").read_text(encoding="utf-8"))
    its = rep["iterations"]
    null = rep["efficacy_null"]
    cands = []
    for it in its:
        t1 = it["report"]["tier1_selectivity"]
        h = t1.get("healthy_control_deltas") or {}
        hv = next(iter(h.values()), 0.0) if isinstance(h, dict) else 0.0
        if isinstance(hv, dict):
            hv = next(iter(hv.values()), 0.0)
        cands.append({"pos": it["report"]["candidate"]["target_edit"]["position"],
                      "lfc": t1["target_delta_sketch_b"],
                      "healthy": float(hv) if isinstance(hv, (int, float)) else 0.0})
    meth = [it["report"]["tier2_hazards"]["cpg_methylation_mean_abs_shift"] for it in its]
    over = [v for v in meth if v > 0.12]
    return {"cands": cands, "abs": [abs(x["lfc"]) for x in cands], "null": null["abs_lfc_values"],
            "meth_closest": min(over) if over else max(meth), "meth_limit": 0.12,
            "cutoff": null["abs_lfc_max"], "n_edits": len(its), "n_null": null["n_random_edits"],
            "n_fail_meth": sum(1 for it in its if not it["report"]["tier2_hazards"]["methylation_hazard_passed"])}


# ================================================================================================
# 0. Intro: present the tool before anything else
# ================================================================================================
def scene_intro(c, t, dur, dark):
    helix(c, 420, 540, 1020, 180, 0.6 + t * 0.45, grow=out_cubic(t / 1.8), dot=15, n=28)
    a, dy = rise(t, 1.0, 1.1)
    w = text(c, "gene-edit-", 840, 400, 112, PAPER, align="left", alpha=a, family=SERIF, bold=False, dy=dy, fit=False)
    text(c, "check", 840 + w, 400, 112, MINT, align="left", alpha=a, family=SERIF, bold=False, dy=dy, fit=False)
    a, dy = rise(t, 1.9, 1.0)
    text(c, "Check a gene edit on a computer,", 844, 530, 46, PAPER, align="left", alpha=a, bold=False, dy=dy)
    text(c, "before anyone spends weeks on it in a lab.", 844, 594, 46, PAPER, align="left", alpha=a, bold=False, dy=dy)
    steps = ["find", "predict", "stress-test", "judge"]
    for i, s in enumerate(steps):
        a2, dy2 = rise(t, 3.3 + i * 0.28, 0.7, 24)
        x = 844 + i * 260
        disc(c, x + 22, 736 + dy2, 22, CORAL, a2 * 0.9)
        text(c, str(i + 1), x + 22, 736 + dy2, 26, PAPER, alpha=a2)
        text(c, s, x + 56, 736, 36, MUTE_D, align="left", alpha=a2, bold=False, dy=dy2, fit=False)
    a3, _ = rise(t, 5.0, 0.9)
    text(c, "four steps, run on real single-cell data", 844, 836, 36, MINT, align="left", alpha=a3, bold=False)


# ================================================================================================
# 1. Find where an edit could work
# ================================================================================================
SEQ_1 = "ACGTTGGCATCGGATCCGTAAGGCTTCAGGACCTGGATGCAAGGTCAGTACGGTACCTGGAAGCCTTGGCATCGGAATCCGTAAGG"


def _find_a(c, t):
    """The locus: a gene body on a ruler, then a window zooms in on the start of the gene."""
    y = 520
    line(c, 200, y, 1720, y, 4, MUTE_L, fade(t, 0.1) * 0.7)
    for i in range(9):
        x = 200 + i * 190
        a = fade(t, 0.15 + i * 0.05, 0.4)
        line(c, x, y - 12, x, y + 12, 3, MUTE_L, a * 0.6)
    a = fade(t, 0.5)
    rrect(c, 760, y - 34, 420, 68, 12, TEAL, a * 0.85)
    text(c, "PDCD1", 970, y, 40, PAPER, alpha=a)
    text(c, "chr2", 200, y + 80, 32, MUTE_L, align="left", alpha=fade(t, 0.3), bold=False)
    text(c, "131 kb around the gene", 1720, y + 80, 32, MUTE_L, align="right", alpha=fade(t, 0.3), bold=False)
    p = out_cubic((t - 1.3) / 1.1)
    if p > 0:
        bx, bw = lerp(200, 880, p), lerp(1520, 180, p)
        for dx in (0, bw):
            line(c, bx + dx, y - 90, bx + dx, y + 90, 4, CORAL, p)
        line(c, bx, y - 90, bx + bw, y - 90, 4, CORAL, p)
        line(c, bx, y + 90, bx + bw, y + 90, 4, CORAL, p)
        text(c, "where the gene starts", bx + bw / 2, y - 130, 34, CORAL, alpha=p, bold=False)


def _find_b(c, t):
    """The model's sensitivity across that window, drawn as a curve over the sequence."""
    x0, x1, base = 200, 1720, 760
    strand(c, x0, x1, base + 60, SEQ_1, False, size=26, alpha=fade(t, 0.1) * 0.8, reveal=out_cubic(t / 1.2))
    xs = np.linspace(0, 1, 260)
    bumps = [(0.50, 0.95, 0.045), (0.43, 0.55, 0.03), (0.60, 0.62, 0.035), (0.29, 0.35, 0.03),
             (0.74, 0.40, 0.04), (0.14, 0.25, 0.03), (0.88, 0.20, 0.03)]
    ys = np.full(len(xs), 0.05)
    for mu, h, s in bumps:
        ys = ys + h * np.exp(-((xs - mu) ** 2) / (2 * s * s))
    ys = ys / ys.max()
    pts = [(x0 + u * (x1 - x0), base - v * 390) for u, v in zip(xs, ys)]
    prog = out_cubic((t - 0.6) / 2.0)
    c.save()
    c.rectangle(0, 0, x0 + prog * (x1 - x0), H)
    c.clip()
    c.move_to(pts[0][0], base)
    for p in pts:
        c.line_to(*p)
    c.line_to(pts[-1][0], base)
    c.close_path()
    c.set_source_rgba(*TEAL, 0.3)
    c.fill()
    c.move_to(*pts[0])
    for p in pts[1:]:
        c.line_to(*p)
    c.set_source_rgba(*TEAL, 1)
    c.set_line_width(5)
    c.set_line_join(cairo.LINE_JOIN_ROUND)
    c.stroke()
    c.restore()
    a = fade(t, 2.3)
    text(c, "the model is most sensitive here", x0 + 0.50 * (x1 - x0), base - 450, 34, INK, alpha=a, bold=False)


def _find_c(c, t):
    """Scan for the places a base editor can actually reach, and pin each one."""
    x0, x1, y = 230, 1690, 560
    step = strand(c, x0, x1, y, SEQ_1, False, size=26, alpha=0.6)
    scan = (t - 0.2) / 2.0
    sites = [i for i in range(10, len(SEQ_1) - 2) if SEQ_1[i] == "G" and SEQ_1[i + 1] == "G"][::2][:6]
    for k, i in enumerate(sites):
        u = (i + 1) / len(SEQ_1)
        if scan < u:
            continue
        p = out_back(min((scan - u) * 4, 1.0))
        rrect(c, x0 + (i - 7) * step - 8, y - 40, step * 9, 80, 10, TEAL, 0.75 * p, stroke=4)
        if k == 0:
            text(c, "PAM", x0 + (i + 0.5) * step, y + 72, 26, TEAL, alpha=p * 0.9, bold=False, fit=False)
            text(c, "where the editor can act", x0 + (i + 1) * step, y - 96, 26, TEAL, alpha=p * 0.9, bold=False, fit=False)
        px = x0 + (i - 4) * step
        line(c, px, y - 52, px, y - 52 - 80 * p, 3, MUTE_L, p * 0.55)
        disc(c, px, y - 142 * p - 2, 14, CORAL, p)
    if 0 <= scan <= 1:
        sx = x0 + scan * (x1 - x0)
        line(c, sx, y - 200, sx, y + 110, 4, INK, 0.25)
    a = fade(t, 1.5)
    text(c, "every spot a base editor can reach", 230, 830, 36, MUTE_L, align="left", alpha=a, bold=False)
    v = counter(t, 1.5, 1.3, 60)
    text(c, f"{v:.0f}", 1690, 810, 96, INK, align="right", alpha=a, family=SERIF, bold=False)
    text(c, "candidates, ranked", 1690, 890, 34, MUTE_L, align="right", alpha=a, bold=False)


def scene_find(c, t, dur, dark):
    header(c, t, dur, 1, "Find where an edit could work", dark)
    cap, a = stage(c, t, [(3.0, _find_a, "one gene, and the stretch of DNA that controls it"),
                          (3.8, _find_b, "input gradients show where the model is sensitive"),
                          (3.9, _find_c, "a base editor needs a PAM site next door")])
    caption(c, cap, a, dark, y=1010)


# ================================================================================================
# 2. Predict what each edit does
# ================================================================================================
def _editor_blob(c, x, y, s, alpha, t):
    """A loose protein shape, so the base editor reads as a molecule rather than a box."""
    c.save()
    c.translate(x, y)
    c.set_source_rgba(*LILAC, alpha * 0.85)
    c.new_path()
    n = 14
    pts = []
    for i in range(n):
        ang = i / n * 2 * math.pi
        rr = s * (1 + 0.16 * math.sin(i * 2.1 + t * 1.4) + 0.08 * math.cos(i * 3.3))
        pts.append((rr * math.cos(ang), rr * math.sin(ang) * 0.82))
    c.move_to(*pts[0])
    for i in range(n):
        p0, p1 = pts[i], pts[(i + 1) % n]
        c.curve_to(p0[0] * 1.08, p0[1] * 1.08, p1[0] * 1.08, p1[1] * 1.08, p1[0], p1[1])
    c.close_path()
    c.fill()
    c.restore()


def _pred_a(c, t):
    """The base editor docks, and one C becomes a T."""
    flipped = t > 1.65
    seq = "GCTTCAGGA" + ("T" if flipped else "C") + "CTGGATGCAAGG"
    x0, x1, y = 420, 1500, 560
    step = strand(c, x0, x1, y, seq, True, size=52, alpha=1.0, highlight={9}, hl_color=CORAL)
    lo, hi = 4, 15
    a = fade(t, 0.2)
    rrect(c, x0 + (lo - 0.6) * step, y - 62, step * (hi - lo + 1.2), 124, 14, MINT, a * 0.5, stroke=4)
    text(c, "guide RNA", x0 + (lo + hi) / 2 * step, y + 112, 32, MINT, alpha=a, bold=False)
    text(c, "PAM", x0 + 19.5 * step, y + 112, 32, GOLD, alpha=fade(t, 0.4), bold=False)
    dock = out_cubic((t - 0.6) / 1.0)
    if dock > 0:
        _editor_blob(c, x0 + 9 * step, lerp(150, 380, dock), 118 * dock, dock, t)
        text(c, "base editor", x0 + 9 * step, lerp(150, 380, dock), 30, FOREST, alpha=dock)
    if flipped:
        pulse = clamp((t - 1.65) / 0.7)
        ring(c, x0 + 9 * step, y, 30 + 100 * pulse, CORAL, 4, 1 - pulse)
        text(c, "C  →  T", x0 + 9 * step, y + 220, 46, CORAL, alpha=fade(t, 1.75), bold=False)
    text(c, "nothing is cut, and nothing is inserted", W / 2, 880, 36, MUTE_D, alpha=fade(t, 2.4), bold=False)


def _pred_b(c, t):
    """The same real edit, read out in two cell types, on one honest axis."""
    cd = max(DATA["cands"], key=lambda q: abs(q["lfc"]) - abs(q["healthy"]))
    pct_t = (2 ** cd["lfc"] - 1) * 100
    pct_h = (2 ** cd["healthy"] - 1) * 100
    text(c, f"chr2:{cd['pos']:,}  ·  the same edit, in two kinds of cell", W / 2, 290, 36, MUTE_D,
         alpha=fade(t, 0.05), bold=False)
    for k, (cx, col, lab, sub) in enumerate(((620, CORAL, "exhausted T cells", "inside the tumour"),
                                             (1300, MINT, "healthy T cells", "normal tissue"))):
        a = fade(t, 0.25 + k * 0.2)
        cell(c, cx, 430, 118, col, alpha=a, fill=1.0, seed=5 + k, t=t, n=20)
        text(c, lab, cx, 582, 38, PAPER, alpha=a)
        text(c, sub, cx, 626, 28, MUTE_D, alpha=a, bold=False)
    # a shared axis, zoomed so a small real change is still readable
    ax0, ax1, ay = 380, 1540, 858
    lo, hi = -15.0, 5.0
    line(c, ax0, ay, ax1, ay, 4, MUTE_D, fade(t, 0.6) * 0.7)
    for v in (-15, -10, -5, 0, 5):
        x = ax0 + (v - lo) / (hi - lo) * (ax1 - ax0)
        line(c, x, ay - 10, x, ay + 10, 3, MUTE_D, fade(t, 0.6) * 0.6)
        text(c, f"{v:+d}%", x, ay + 48, 26, MUTE_D, alpha=fade(t, 0.6) * 0.8, bold=False, fit=False)
    zero = ax0 + (0 - lo) / (hi - lo) * (ax1 - ax0)
    line(c, zero, ay - 175, zero, ay + 24, 3, MUTE_D, fade(t, 0.6) * 0.5, dash=[10, 8])
    prog = out_cubic((t - 1.2) / 1.6)
    for k, (pct, col, yy) in enumerate(((pct_t, CORAL, ay - 148), (pct_h, MINT, ay - 86))):
        a = fade(t, 1.1 + k * 0.15)
        v = pct * prog
        x = ax0 + (v - lo) / (hi - lo) * (ax1 - ax0)
        rrect(c, min(x, zero), yy - 15, abs(x - zero), 30, 15, col, a * 0.85)
        disc(c, x, yy, 17, col, a)
        text(c, f"{v:+.1f}%", x - 28, yy, 32, col, align="right", alpha=a, bold=False, fit=False)
    text(c, "PD-1 the cell makes, against no edit at all", W / 2, 988, 32, MUTE_D, alpha=fade(t, 2.2), bold=False)


def _pred_c(c, t):
    """All the verified candidates on one scale: the gap between the two dots is the selectivity."""
    rows = DATA["cands"][:5]
    text(c, "every candidate, scored the same way", W / 2, 285, 36, MUTE_D, alpha=fade(t, 0.05), bold=False)
    ax0, ax1 = 620, 1620
    lo, hi = -32.0, 8.0

    def px(v):
        return ax0 + (v - lo) / (hi - lo) * (ax1 - ax0)

    zero = px(0)
    aa = fade(t, 0.2)
    line(c, zero, 330, zero, 880, 3, MUTE_D, aa * 0.45, dash=[10, 8])
    for v in (-30, -20, -10, 0):
        text(c, f"{v:+d}%", px(v), 912, 26, MUTE_D, alpha=aa * 0.7, bold=False, fit=False)
    for i, cd in enumerate(rows):
        y = 400 + i * 100
        a, dy = rise(t, 0.3 + i * 0.2, 0.7, 26)
        text(c, f"chr2:{cd['pos']:,}", 230, y + dy, 30, PAPER, align="left", alpha=a, family=MONO, bold=False, fit=False)
        prog = out_cubic((t - 0.55 - i * 0.2) / 0.9)
        xt, xh = px((2 ** cd["lfc"] - 1) * 100 * prog), px((2 ** cd["healthy"] - 1) * 100 * prog)
        line(c, xt, y + dy, xh, y + dy, 4, PAPER, a * 0.28)
        disc(c, xh, y + dy, 14, MINT, a)
        disc(c, xt, y + dy, 16, CORAL, a)
        gap = abs(cd["lfc"]) - abs(cd["healthy"])
        a2 = fade(t, 1.7 + i * 0.16)
        text(c, "selective" if gap > 0.02 else "healthy cells move too", 1800, y + dy, 28,
             MINT if gap > 0.02 else MUTE_D, align="right", alpha=a2, bold=False, fit=False)
    a3 = fade(t, 0.9)
    disc(c, 640, 980, 14, CORAL, a3)
    text(c, "exhausted T cells", 668, 980, 30, CORAL, align="left", alpha=a3, bold=False, fit=False)
    disc(c, 1000, 980, 12, MINT, a3)
    text(c, "healthy T cells", 1028, 980, 30, MINT, align="left", alpha=a3, bold=False, fit=False)


def scene_predict(c, t, dur, dark):
    header(c, t, dur, 2, "Predict what each edit does", dark)
    cap, a = stage(c, t, [(4.3, _pred_a, "a single base, swapped in place"),
                          (4.4, _pred_b, "AlphaGenome 2 reads it in each cell type"),
                          (4.2, _pred_c, "selective means: moves the target, leaves the rest")])
    caption(c, cap, a, dark, y=1050)


# ================================================================================================
# 3. Stress-test every edit  (shown, not listed)
# ================================================================================================
def _stress_splice(c, t):
    """Exons joined by an arc, until a cryptic site drags the arc to the wrong place."""
    y = 560
    exons = [(360, 240), (760, 180), (1180, 300)]
    for i, (x, w) in enumerate(exons):
        a, dy = rise(t, 0.1 + i * 0.12, 0.6, 24)
        rrect(c, x, y - 36 + dy, w, 72, 14, TEAL, a * 0.9)
    line(c, 600, y, 760, y, 4, MUTE_L, fade(t, 0.3) * 0.5)
    line(c, 940, y, 1180, y, 4, MUTE_L, fade(t, 0.3) * 0.5)
    bad = clamp((t - 1.7) / 0.9)
    for k, (sx, ex, col, al) in enumerate(((600, 760, INK, 1 - bad), (600, 1060, CORAL, bad))):
        if al <= 0.01:
            continue
        c.set_source_rgba(*col, al * 0.9)
        c.set_line_width(5)
        c.move_to(sx, y - 36)
        c.curve_to(sx + (ex - sx) * 0.3, y - 190, ex - (ex - sx) * 0.3, y - 190, ex, y - 36)
        c.stroke()
    if bad > 0:
        disc(c, 1060, y, 14, CORAL, bad)
        text(c, "a new splice site appears", 1060, y + 150, 36, CORAL, alpha=bad, bold=False)
    text(c, "splicing", 360, y - 150, 54, INK, align="left", alpha=fade(t, 0.1), family=SERIF, bold=False)


def _stress_methyl(c, t):
    """CpG tags blooming along the strand, with a gauge that climbs past the fixed limit."""
    y = 740
    line(c, 260, y, 1300, y, 6, INK, fade(t, 0.05) * 0.6)
    for i in range(18):
        x = 290 + i * 58
        p = out_back(clamp((t - (0.4 + abs(i - 9) * 0.09)) / 0.5))
        if p <= 0:
            continue
        bob = 4 * math.sin(t * 2 + i)
        line(c, x, y, x, y - 42 + bob, 3, INK, p * 0.5)
        disc(c, x, y - 56 + bob, 14 * p, CORAL, p * 0.9)
    text(c, "methylation", 260, 400, 54, INK, align="left", alpha=fade(t, 0.05), family=SERIF, bold=False)
    text(c, "chemical tags that help switch genes on and off", 260, 470, 32, MUTE_L, align="left",
         alpha=fade(t, 0.3), bold=False)
    # the gauge: 0 at the bottom, the limit as a dashed line, the real shift climbing past it
    gx, gy0, gy1 = 1480, y, 330
    vmax = 0.30
    a = fade(t, 0.8)
    rrect(c, gx, gy1, 130, gy0 - gy1, 18, INK, a * 0.08)
    ylim = gy0 - (DATA["meth_limit"] / vmax) * (gy0 - gy1)
    v = DATA["meth_closest"] * out_cubic((t - 1.1) / 1.6)
    yv = gy0 - (v / vmax) * (gy0 - gy1)
    over = v > DATA["meth_limit"]
    rrect(c, gx, yv, 130, gy0 - yv, 18, CORAL if over else TEAL, a * 0.9)
    line(c, gx - 50, ylim, gx + 180, ylim, 4, INK, a * 0.8, dash=[12, 10])
    text(c, f"limit  {DATA['meth_limit']:.2f}", gx - 60, ylim, 30, INK, align="right", alpha=a, bold=False, fit=False)
    if over:
        ao = clamp((v - DATA["meth_limit"]) / 0.02)
        text(c, f"{v:.3f}", gx + 65, yv - 46, 44, CORAL, alpha=ao, family=SERIF, bold=False, fit=False)
        text(c, "the closest any candidate got", gx + 65, 270, 28, CORAL, alpha=ao, bold=False, fit=False)


def _stress_rest(c, t):
    """The remaining four hazards as quick, wordless vignettes."""
    cells = [(0, "neighbouring genes"), (1, "bystander edits"), (2, "other tissues"), (3, "healthy cells")]
    for i, (k, lab) in enumerate(cells):
        x = 300 + i * 440
        a, dy = rise(t, 0.15 + i * 0.3, 0.7, 30)
        if a <= 0:
            continue
        cy = 540 + dy
        if k == 0:
            for j, dx in enumerate((-110, 0, 110)):
                rrect(c, x + dx - 42, cy - 26, 84, 56, 12, TEAL if j == 1 else MUTE_L, a * (0.95 if j == 1 else 0.3))
            r = ((t - 0.15 - i * 0.3) % 1.8) / 1.8
            ring(c, x, cy, 30 + 90 * r, TEAL, 3, a * (1 - r) * 0.8)
        elif k == 1:
            s = "ATCGGCA"
            for j, ch in enumerate(s):
                hit = j in (3, 4)
                text(c, ch, x - 108 + j * 36, cy, 40, CORAL if hit else MUTE_L, alpha=a * (1 if hit else 0.6),
                     family=MONO, fit=False)
            rrect(c, x - 20, cy - 36, 108, 72, 10, CORAL, a * 0.7, stroke=3)
        elif k == 2:
            for j, cc in enumerate((LILAC, (0.40, 0.65, 0.80), (0.92, 0.55, 0.52), (0.55, 0.78, 0.55), GOLD)):
                disc(c, x - 104 + j * 52, cy + 6 * math.sin(t * 1.6 + j), 22, cc, a * 0.9)
        else:
            for dx in (-62, 62):
                ring(c, x + dx, cy, 44, TEAL, 5, a)
                disc(c, x + dx - 10, cy - 8, 15, TEAL, a * 0.35)
        text(c, lab, x, cy + 170, 34, MUTE_L, alpha=a, bold=False)
        res = clamp((t - 1.4 - i * 0.22) / 0.5)
        tick(c, x, cy + 250, out_back(res), TEAL)
    text(c, "and four more", 300, 300, 54, INK, align="left", alpha=fade(t, 0.05), family=SERIF, bold=False)


def _stress_sum(c, t):
    """All six at once, with the one that actually fails called out."""
    names = ["splicing", "methylation", "neighbours", "bystanders", "tissues", "healthy cells"]
    pw, gx = 250, 26
    x0 = (W - (6 * pw + 5 * gx)) / 2
    for i, n in enumerate(names):
        x = x0 + i * (pw + gx)
        a, dy = rise(t, 0.1 + i * 0.1, 0.6, 26)
        bad = i == 1
        rrect(c, x, 440 + dy, pw, 250, 24, CARD, a)
        if bad and t > 1.3:
            rrect(c, x, 440 + dy, pw, 250, 24, CORAL, a * fade(t, 1.3), stroke=5)
        res = clamp((t - 0.8 - i * 0.1) / 0.5)
        if bad:
            cross(c, x + pw / 2, 530 + dy, out_back(clamp((t - 1.3) / 0.5)), CORAL, 10)
        else:
            tick(c, x + pw / 2, 530 + dy, out_back(res), TEAL, 10)
        text(c, n, x + pw / 2, 630 + dy, 30, INK, alpha=a, bold=False)
    a2 = fade(t, 2.0)
    text(c, f"{DATA['n_fail_meth']} of {DATA['n_edits']} candidates died right here", W / 2, 820, 48, INK,
         alpha=a2, family=SERIF, bold=False)


def scene_stress(c, t, dur, dark):
    header(c, t, dur, 3, "Stress-test every edit", dark)
    cap, a = stage(c, t, [(3.4, _stress_splice, "does it break how the gene is read?"),
                          (3.4, _stress_methyl, "does it disturb the chemical tags on the DNA?"),
                          (3.2, _stress_rest, "and the damage an edit does next door"),
                          (3.0, _stress_sum, "six limits, all fixed before any result was seen")])
    caption(c, cap, a, dark, y=940)


# ================================================================================================
# 4. Judge on patients it never saw  (slow)
# ================================================================================================
def _judge_a(c, t):
    """One cohort of patients splits in two, and the halves never meet again."""
    rng = np.random.default_rng(21)
    pts = [(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(46)]
    split = out_cubic((t - 1.6) / 1.8)
    for i, (ux, uy) in enumerate(pts):
        a = fade(t, 0.15 + i * 0.02, 0.5)
        grp = i % 2
        cx = lerp(960, 600 if grp == 0 else 1320, split)
        col = PAPER if split < 0.25 else (MINT if grp == 0 else GOLD)
        disc(c, cx + ux * lerp(300, 180, split), 560 + uy * lerp(230, 190, split), 15,
             col, a * 0.9)
    text(c, "62 patients", 960, 250, 44, PAPER, alpha=fade(t, 0.4) * (1 - split), bold=False)
    for k, (cx, col, lab) in enumerate(((600, MINT, "group A"), (1320, GOLD, "group B"))):
        a = split
        ring(c, cx, 560, 230, col, 3, a * 0.35)
        text(c, lab, cx, 860, 42, col, alpha=a)
    line(c, 960, 320, 960, 800, 4, CORAL, split * 0.8, dash=[14, 12])
    text(c, "never mixed", 960, 280, 34, CORAL, alpha=split, bold=False)


def _judge_b(c, t):
    """An edit is proposed on one half and re-measured on the other."""
    for k, (cx, col, lab) in enumerate(((560, MINT, "group A  ·  proposes"), (1360, GOLD, "group B  ·  judges"))):
        ring(c, cx, 560, 210, col, 3, fade(t, 0.1) * 0.5)
        rng = np.random.default_rng(31 + k)
        for j in range(18):
            ang, rr = rng.uniform(0, 2 * math.pi), 180 * math.sqrt(rng.uniform(0.05, 1))
            disc(c, cx + rr * math.cos(ang), 560 + rr * math.sin(ang), 14, col, fade(t, 0.1) * 0.75)
        text(c, lab, cx, 830, 38, PAPER, alpha=fade(t, 0.2), bold=False)
    found = out_back(clamp((t - 0.9) / 0.7))
    if found > 0:
        disc(c, 560, 560, 26 * found, CORAL, 1)
        text(c, "best edit found here", 560, 300, 32, CORAL, alpha=clamp(found), bold=False)
    travel = out_cubic((t - 1.9) / 1.4)
    if travel > 0:
        x = lerp(560, 1360, travel)
        line(c, 560, 560, x, 560, 4, CORAL, 0.5)
        disc(c, x, 560, 26, CORAL, 1)
    if t > 3.3:
        a = fade(t, 3.3)
        text(c, "measured again, from scratch", 1360, 300, 32, GOLD, alpha=a, bold=False)
        rrect(c, 1160, 918, 400, 26, 13, PAPER, a * 0.12)
        rrect(c, 1160, 918, 400 * (2 ** DATA["cands"][0]["lfc"]), 26, 13, CORAL, a)
        text(c, "the effect holds up", 1360, 982, 34, MINT, alpha=fade(t, 4.0), bold=False)


def _judge_c(c, t):
    """And it still has to beat a crowd of random edits."""
    ax0, ax1, ay = 240, 1680, 720
    line(c, ax0, ay, ax1, ay, 4, MUTE_D, fade(t, 0.2) * 0.8)
    xmax = max(DATA["abs"]) * 1.12

    def px(v):
        return ax0 + 40 + math.sqrt(min(v, xmax) / xmax) * (ax1 - ax0 - 80)

    jit = np.random.default_rng(77).normal(0, 11, len(DATA["null"]))
    for i, v in enumerate(DATA["null"]):
        a, dy = rise(t, 0.5 + i * 0.045, 0.55, 36)
        disc(c, px(v) + jit[i], ay - 34 - 30 * (i % 5) + dy, 13, MINT, a * 0.75)
    text(c, f"{DATA['n_null']} random edits", ax0, 330, 40, MINT, align="left", alpha=fade(t, 0.6), bold=False)
    a = fade(t, 2.4)
    cx = px(DATA["cutoff"])
    line(c, cx, ay + 34, cx, ay - 300, 4, GOLD, a)
    text(c, "the bar: the best of them", cx + 18, ay - 330, 34, GOLD, align="left", alpha=a, bold=False)
    for i, v in enumerate(DATA["abs"]):
        a2, dy = rise(t, 3.1 + i * 0.16, 0.7, 36)
        over = v > DATA["cutoff"]
        disc(c, px(v), ay - 44 - 42 * (i % 3) + dy, 21, CORAL if over else MUTE_D, a2 * (1 if over else 0.55))
    text(c, "only what clears it counts", ax1, ay + 120, 36, PAPER, align="right", alpha=fade(t, 4.4), bold=False)


def scene_judge(c, t, dur, dark):
    header(c, t, dur, 4, "Judge on patients it never saw", dark)
    cap, a = stage(c, t, [(4.6, _judge_a, "the cohort is split by patient, once, before anything runs"),
                          (5.4, _judge_b, "proposed on one half, re-measured on the other"),
                          (5.0, _judge_c, "and it still has to beat chance")])
    caption(c, cap, a, dark, y=1030)


# ================================================================================================
# 5. What it found
# ================================================================================================
def scene_stats(c, t, dur, dark):
    a, dy = rise(t, 0.15, 0.8)
    text(c, "What the first real run found", 150, 150, 70, INK, align="left", alpha=a, family=SERIF, bold=False, dy=dy)
    tiles = [(150, 430, 62, "patients", 0.6, INK), (1030, 430, 12738, "exhausted T cells", 1.0, INK),
             (150, 740, DATA["n_edits"], "edits fully verified", 1.4, INK), (1030, 740, 0, "passed every check", 1.8, CORAL)]
    for x, y, v, lab, t0, col in tiles:
        a, dy = rise(t, t0, 0.9, 28)
        text(c, f"{v * out_cubic((t - t0) / 1.1):,.0f}", x, y, 150, col, align="left", alpha=a,
             family=SERIF, bold=False, dy=dy, fit=False)
        text(c, lab, x + 4, y + 100, 38, MUTE_L, align="left", alpha=a, bold=False, dy=dy, fit=False)
    a, dy = rise(t, 2.6, 0.9)
    text(c, f"{DATA['n_fail_meth']} of {DATA['n_edits']} were stopped by methylation. A negative result, stated plainly.",
         150, 960, 44, INK, align="left", alpha=a, family=SERIF, bold=False, dy=dy)


# ================================================================================================
# 6. Use it on your own data  (the detailed part)
# ================================================================================================
def _term(c, x, y, w, h, lines, shown, *, alpha=1.0, title=None):
    rrect(c, x, y, w, h, 24, TERM, alpha)
    rrect(c, x, y, w, h, 24, MINT, alpha * 0.18, stroke=2)
    for k, col in enumerate((CORAL, GOLD, MINT)):
        disc(c, x + 42 + k * 34, y + 42, 9, col, alpha * 0.9)
    if title:
        text(c, title, x + w - 40, y + 42, 26, MUTE_D, align="right", alpha=alpha * 0.8, bold=False, fit=False)
    yy = y + 110
    left = shown
    for s in lines:
        n = min(len(s), left)
        left -= max(n, 0)
        if n <= 0:
            break
        indent = len(s) - len(s.lstrip())
        text(c, s[:n].lstrip(), x + 44 + indent * 20, yy, 34,
             MINT if s.startswith("$") else PAPER, align="left", alpha=alpha, family=MONO, bold=False, fit=False)
        yy += 62
    return yy


def _use_a(c, t):
    """Install it, then let it check its own setup."""
    lines = ["$ pip install gene-edit-check", "$ gene-edit-check doctor"]
    total = sum(len(s) for s in lines)
    _term(c, 180, 300, 880, 320, lines, int(total * clamp((t - 0.2) / 1.6)), title="any terminal")
    rows = [("API key", "found"), ("reference conditions", "2,799 shipped"), ("cell atlas", "downloads on first run"),
            ("cache", "per-user folder")]
    for i, (k, v) in enumerate(rows):
        a, dy = rise(t, 1.9 + i * 0.22, 0.6, 24)
        y = 340 + i * 90 + dy
        rrect(c, 1140, y - 34, 620, 72, 16, PAPER, a * 0.07)
        tick(c, 1190, y, out_back(clamp((t - 2.0 - i * 0.22) / 0.5)), MINT, 7)
        text(c, k, 1240, y, 32, PAPER, align="left", alpha=a, bold=False, fit=False)
        text(c, v, 1740, y, 28, MUTE_D, align="right", alpha=a, bold=False, fit=False)
    text(c, "it tells you what is missing before you spend an API call", W / 2, 760, 34, MUTE_D,
         alpha=fade(t, 3.0), bold=False)


def _use_b(c, t):
    """Point it at your own single-cell file, and say which cells are which."""
    lines = ["$ gene-edit-check run-h5ad \\", "    --h5ad my_cells.h5ad \\", "    --gene-symbol PDCD1 \\",
             "    --target-query \"state == 'exhausted'\" \\", "    --healthy-query \"state == 'naive'\" \\",
             "    --scan-region 14300:15300"]
    total = sum(len(s) for s in lines)
    _term(c, 180, 250, 1080, 520, lines, int(total * clamp((t - 0.15) / 2.6)), title="your data")
    notes = [("your own cells", 0.0), ("any gene", 1.0), ("the cells to change", 2.0),
             ("the cells to protect", 3.0), ("or aim at an enhancer", 4.0)]
    for i, (s, row) in enumerate(notes):
        a, dx = rise(t, 1.5 + i * 0.3, 0.6, -30)
        y = 360 + (row + 1) * 62
        line(c, 1270, y, 1330, y, 3, CORAL, a * 0.6)
        text(c, s, 1350 + dx, y, 32, CORAL, align="left", alpha=a, bold=False, fit=False)
    text(c, "no code to write, and nothing hard-coded to PD-1", W / 2, 880, 34, MUTE_D, alpha=fade(t, 3.2), bold=False)


def _use_c(c, t):
    """What comes back out."""
    text(c, "what you get back", W / 2, 290, 44, PAPER, alpha=fade(t, 0.05), family=SERIF, bold=False)
    cards = [
        (200, "report.json", ["every candidate", "every check", "the exact thresholds", "model version + seed"]),
        (740, "figures", None),
        (1280, "candidate", ["chr2:241,859,845", "A → G  ·  ABE8e", "protospacer + PAM", "Atlas deep link"]),
    ]
    for i, (x, title_, rows) in enumerate(cards):
        a, dy = rise(t, 0.3 + i * 0.35, 0.8, 40)
        rrect(c, x, 360 + dy, 440, 440, 26, PAPER, a * 0.07)
        rrect(c, x, 360 + dy, 440, 440, 26, MINT, a * 0.22, stroke=2)
        text(c, title_, x + 36, 420 + dy, 36, MINT, align="left", alpha=a, family=MONO, bold=False, fit=False)
        if rows:
            for j, r in enumerate(rows):
                a2 = fade(t, 0.7 + i * 0.35 + j * 0.12)
                disc(c, x + 48, 500 + j * 62 + dy, 5, MUTE_D, a2)
                text(c, r, x + 72, 500 + j * 62 + dy, 30, PAPER, align="left", alpha=a2, bold=False, fit=False)
        else:
            gx, gy = x + 44, 490 + dy
            for j in range(4):
                for k in range(6):
                    a2 = fade(t, 0.8 + j * 0.08 + k * 0.03, 0.4)
                    good = not (k == 3 and j != 1)
                    if good:
                        tick(c, gx + 30 + k * 58, gy + j * 56, out_back(a2) * 0.7, TEAL, 6)
                    else:
                        cross(c, gx + 30 + k * 58, gy + j * 56, out_back(a2) * 0.7, CORAL, 6)
            text(c, "light and dark, print-ready", x + 36, 740 + dy, 26, MUTE_D, align="left",
                 alpha=fade(t, 1.6), bold=False, fit=False)
    text(c, "and a one-page figure you can put straight in a slide deck", W / 2, 880, 34, MUTE_D,
         alpha=fade(t, 2.2), bold=False)


def _use_d(c, t):
    """The other ways in, and where it runs."""
    text(c, "three ways to drive it", W / 2, 290, 44, PAPER, alpha=fade(t, 0.05), family=SERIF, bold=False)
    ways = [("command line", ["gene-edit-check optimize", "gene-edit-check sweep-coupling"]),
            ("Python", ["from gene_edit_check import", "    run_scenario"]),
            ("AI agents", ["MCP server + agent skill", "Claude, Cursor, anything"])]
    for i, (title_, rows) in enumerate(ways):
        a, dy = rise(t, 0.3 + i * 0.3, 0.8, 36)
        x = 170 + i * 540
        rrect(c, x, 360 + dy, 480, 260, 26, PAPER, a * 0.07)
        text(c, title_, x + 36, 420 + dy, 38, MINT, align="left", alpha=a, bold=False, fit=False)
        for j, r in enumerate(rows):
            text(c, r, x + 36, 500 + j * 50 + dy, 26, MUTE_D, align="left",
                 alpha=fade(t, 0.7 + i * 0.3 + j * 0.12), family=MONO, bold=False, fit=False)
    chips = ["macOS", "Linux", "Windows", "Python 3.11+", "MIT licensed"]
    for i, s in enumerate(chips):
        a, dy = rise(t, 1.7 + i * 0.16, 0.6, 24)
        x = 230 + i * 300
        rrect(c, x, 740 + dy, 260, 86, 43, MINT, a * 0.6, stroke=2)
        text(c, s, x + 130, 783 + dy, 30, PAPER, alpha=a, fit=False)
    text(c, "every run writes down its seed, its thresholds and the model version", W / 2, 920, 34, MUTE_D,
         alpha=fade(t, 2.6), bold=False)


def scene_use(c, t, dur, dark):
    header(c, t, dur, 0, "Use it on your own data", dark)
    cap, a = stage(c, t, [(4.2, _use_a, "install it, and check the setup"),
                          (4.6, _use_b, "point it at your cells and name the two states"),
                          (4.2, _use_c, "it writes the evidence down, not just a verdict"),
                          (4.1, _use_d, "a terminal, a Python import, or an agent")])
    caption(c, cap, a, dark, y=1020)


def scene_end(c, t, dur, dark):
    helix(c, 460, 540, 900, 160, 0.8 + t * 0.45, grow=out_cubic(t / 1.4), dot=14, n=24)
    a, dy = rise(t, 0.4, 1.0)
    w = text(c, "gene-edit-", 880, 400, 106, PAPER, align="left", alpha=a, family=SERIF, bold=False, dy=dy, fit=False)
    text(c, "check", 880 + w, 400, 106, MINT, align="left", alpha=a, family=SERIF, bold=False, dy=dy, fit=False)
    a, dy = rise(t, 1.1, 0.9)
    text(c, "Check a gene edit on a computer", 884, 520, 42, PAPER, align="left", alpha=a, bold=False, dy=dy)
    text(c, "before anyone spends weeks on it in a lab.", 884, 582, 42, PAPER, align="left", alpha=a, bold=False, dy=dy)
    s = "pip install gene-edit-check"
    n = int(len(s) * clamp((t - 1.9) / 1.1))
    text(c, s[:n] + ("_" if int(t * 2) % 2 == 0 else ""), 884, 730, 46, GOLD, align="left", family=MONO, bold=False)
    text(c, "github.com/aravmendiratta/gene-edit-check", 884, 820, 34, MUTE_D, align="left",
         alpha=fade(t, 2.9), bold=False)


SCENE_FN = {"intro": scene_intro, "find": scene_find, "predict": scene_predict, "stress": scene_stress,
            "judge": scene_judge, "stats": scene_stats, "use": scene_use, "end": scene_end}


def draw_scene(c, idx, t):
    name, bars, dark = SCENES[idx]
    dur = bars * BAR
    c.set_source_rgb(*(FOREST if dark else PAPER))
    c.paint()
    c.save()
    z = 1.0 + 0.016 * clamp(t / dur)  # a very slow push-in, just enough to stop the frame feeling static
    c.translate(W / 2, H / 2)
    c.scale(z, z)
    c.translate(-W / 2, -H / 2)
    SCENE_FN[name](c, t, dur, dark)
    c.restore()


def render(t, surf):
    c = cairo.Context(surf)
    c.scale(OW / W, OH / H)
    idx = max(i for i in range(len(SCENES)) if STARTS[i] <= t + 1e-9)
    local = t - STARTS[idx]
    if idx > 0 and local < 0.55:
        draw_scene(c, idx - 1, SCENES[idx - 1][1] * BAR + local)
        c.push_group()
        draw_scene(c, idx, local)
        c.pop_group_to_source()
        c.paint_with_alpha(in_out(local / 0.55))
    else:
        draw_scene(c, idx, local)


# ---- music ---------------------------------------------------------------------------------------
SR = 44100


def make_music() -> np.ndarray:
    """112 BPM, A minor. Pads and bells carry it; a soft kick, shaker and arpeggio keep it moving."""
    from scipy.signal import butter, fftconvolve, lfilter

    n = int(TOTAL * SR) + SR
    dry = np.zeros((n, 2))
    wet = np.zeros((n, 2))
    rng = np.random.default_rng(12)

    def tt(d):
        return np.arange(int(d * SR)) / SR

    def filt(x, kind, fc, order=2):
        b, a = butter(order, fc / (SR / 2), btype=kind)
        return lfilter(b, a, x)

    def put(buf, t0, sig, gl=1.0, gr=None):
        i = int(t0 * SR)
        if i < 0 or i >= n:
            return
        m = min(len(sig), n - i)
        buf[i : i + m, 0] += sig[:m] * gl
        buf[i : i + m, 1] += sig[:m] * (gl if gr is None else gr)

    def saw(f, d, de=0.0):
        return 2 * ((tt(d) * f * (1 + de)) % 1.0) - 1

    def bell(f, d=1.8):
        t = tt(d)
        s = (np.sin(2 * np.pi * f * t) + 0.3 * np.sin(2 * np.pi * 2 * f * t) * np.exp(-t * 5)
             + 0.12 * np.sin(2 * np.pi * 3.01 * f * t) * np.exp(-t * 8))
        return s * np.exp(-t * 2.8) * np.minimum(1, t * 250) * 0.26

    def pluck(f, d=0.5):
        t = tt(d)
        return filt(saw(f, d) * 0.5 + np.sin(2 * np.pi * f * t) * 0.5, "low", 2600) * np.exp(-t * 7) * 0.11

    def kick(d=0.34):
        t = tt(d)
        return np.sin(2 * np.pi * np.cumsum(50 + 70 * np.exp(-t * 28)) / SR) * np.exp(-t * 10) * 0.45

    def shaker(d=0.07):
        t = tt(d)
        return filt(rng.normal(0, 1, len(t)), "high", 6800) * np.exp(-t * 55) * 0.05

    chords = [(110.0, [220.0, 261.63, 329.63]), (87.31, [174.61, 220.0, 261.63]),
              (130.81, [196.0, 261.63, 329.63]), (98.0, [196.0, 246.94, 293.66])]
    bars = int(round(TOTAL / BAR)) + 1
    scene_bars = {int(round(s / BAR)) for s in STARTS}

    for bar in range(bars):
        root, triad = chords[bar % 4]
        t_bar = bar * BAR
        # pad
        d = BAR + 1.2
        env = np.minimum(1, tt(d) / 0.8) * np.minimum(1, (d - tt(d)) / 1.0)
        sl = sum(saw(f, d, de) for f in triad + [triad[0] * 2] for de in (-0.003, 0.002)) / 8
        sr_ = sum(saw(f, d, de) for f in triad + [triad[0] * 2] for de in (0.003, -0.002)) / 8
        put(wet, t_bar, filt(sl, "low", 1200) * env * 0.5, 1.0, 0.0)
        put(wet, t_bar, filt(sr_, "low", 1200) * env * 0.5, 0.0, 1.0)
        # sub bass, from the second bar on
        if bar >= 1:
            put(dry, t_bar, np.sin(2 * np.pi * root * tt(BAR * 1.05)) * np.exp(-tt(BAR * 1.05) * 1.3) * 0.26)
        # bells answering across the field
        for k, b in enumerate([0.0, 1.5, 3.0] if bar % 2 == 0 else [0.5, 2.0, 3.5]):
            f = triad[(bar + k) % 3] * (2 if k % 2 == 0 else 4)
            tb = t_bar + b * BEAT
            pan = 0.35 + 0.3 * (k % 2)
            put(wet, tb, bell(f), 1 - pan, pan)
            put(dry, tb, bell(f) * 0.22, 1 - pan, pan)
        # the pulse: kick on 1 and 3, shaker on the offbeats, a light 16th arpeggio
        if bar >= 2:
            put(dry, t_bar, kick(), 0.9)
            put(dry, t_bar + 2 * BEAT, kick(), 0.65)
            for e in range(8):
                put(dry, t_bar + e * BEAT / 2 + BEAT / 4, shaker(), 0.75, 0.9)
        if bar >= 3:
            for s16 in range(16):
                f = triad[(s16 + bar) % 3] * (2 if s16 % 4 == 1 else 1)
                tp = t_bar + s16 * BEAT / 4
                pan = 0.4 + 0.25 * (s16 % 2)
                put(dry, tp, pluck(f), 1 - pan, pan)
    # a soft swell into every scene change
    for s in range(1, len(SCENES)):
        t0 = float(STARTS[s])
        d = 1.4
        put(wet, t0 - d, filt(rng.normal(0, 1, int(d * SR)), "low", 2400) * (tt(d) / d) ** 2 * 0.09)
        put(dry, t0, np.sin(2 * np.pi * 55 * tt(1.6)) * np.exp(-tt(1.6) * 2.4) * 0.26)

    ir_t = tt(2.4)
    ir = filt(rng.normal(0, 1, len(ir_t)), "low", 3600) * np.exp(-ir_t * 2.6)
    ir /= np.abs(ir).sum() / 13
    rev = np.stack([fftconvolve(wet[:, 0], ir)[:n], fftconvolve(wet[:, 1], np.roll(ir, 197))[:n]], axis=1)
    mix = (dry + rev * 0.9 + wet * 0.3)[: int(TOTAL * SR)]
    m = len(mix)
    mix *= (np.minimum(1, np.arange(m) / SR / 0.8) * np.minimum(1, (m - np.arange(m)) / SR / 2.2))[:, None]
    mix = np.tanh(mix * 0.95)
    mix /= max(1e-6, np.abs(mix).max()) / 0.74
    return (mix * 32767).astype(np.int16)


# ---- encode ----------------------------------------------------------------------------------
def main(out: str = "assets/launch.mp4") -> None:
    import wave

    import imageio_ffmpeg

    DATA.update(load_data())
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    wav = Path(out).with_suffix(".wav")
    with wave.open(str(wav), "wb") as w:
        w.setnchannels(2)
        w.setsampwidth(2)
        w.setframerate(SR)
        w.writeframes(make_music().tobytes())
    cmd = [imageio_ffmpeg.get_ffmpeg_exe(), "-y", "-loglevel", "error",
           "-f", "rawvideo", "-pix_fmt", "bgra", "-s", f"{OW}x{OH}", "-r", str(FPS), "-i", "-", "-i", str(wav),
           "-c:v", "libx264", "-profile:v", "main", "-level", "4.0", "-preset", "slow", "-crf", "22",
           "-maxrate", "3M", "-bufsize", "6M", "-bf", "0", "-g", str(FPS), "-keyint_min", str(FPS),
           "-sc_threshold", "0", "-pix_fmt", "yuv420p", "-movflags", "+faststart",
           "-c:a", "aac", "-b:a", "160k", "-ar", str(SR), "-shortest", out]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    surf = cairo.ImageSurface(cairo.FORMAT_ARGB32, OW, OH)
    total = int(round(TOTAL * FPS))
    for f in range(total):
        render(f / FPS, surf)
        proc.stdin.write(bytes(surf.get_data()))
    proc.stdin.close()
    proc.wait()
    wav.unlink()
    print(f"wrote {out}: {total} frames, {total / FPS:.1f}s, with sound")


def preview(times: list[float], path: str) -> None:
    from PIL import Image

    DATA.update(load_data())
    surf = cairo.ImageSurface(cairo.FORMAT_ARGB32, OW, OH)
    tiles = []
    for t in times:
        render(t, surf)
        buf = np.ndarray((OH, OW, 4), np.uint8, surf.get_data())
        tiles.append(Image.fromarray(np.ascontiguousarray(buf[:, :, [2, 1, 0]])))
    cols = 2
    sheet = Image.new("RGB", (OW * cols, OH * math.ceil(len(tiles) / cols)))
    for i, im in enumerate(tiles):
        sheet.paste(im, ((i % cols) * OW, (i // cols) * OH))
    sheet.save(path)


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--preview":
        preview([float(x) for x in sys.argv[3:]], sys.argv[2])
    else:
        main(*sys.argv[1:2])
