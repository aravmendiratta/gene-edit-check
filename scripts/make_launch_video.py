"""The gene-edit-check launch film: a Neo-Bauhaus explainer, drawn and scored entirely in code.

    pip install pycairo imageio-ffmpeg numpy scipy
    PYTHONPATH=src python scripts/make_launch_video.py            # writes assets/launch.mp4 (1280x720)
    PYTHONPATH=src python scripts/make_launch_video.py --preview out.png 12 34 56

Design: circle, square and rule only. Primary red, blue and yellow on warm paper or near-black, asymmetric
compositions with shapes bleeding off the frame, and a different geometric transition between every chapter.
The flow is title, problem, use case, benefits, outro.

Sound: a 108 BPM bed plus a sound effect on every element that lands, all synthesised here, so nothing needs
licensing. `build_events()` mirrors the animation schedules held in `TIMING`; the two must be edited together,
which is why every schedule is a named constant rather than a literal buried in a draw call.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import cairo
import numpy as np

W, H = 1920, 1080            # layout grid
OW, OH, FPS = 1280, 720, 30  # render size
BPM = 108
BEAT = 60 / BPM
BAR = BEAT * 4
OUT = Path("outputs")

# Bauhaus primaries, knocked back so they read as pigment rather than screen light.
PAPER = (0.945, 0.929, 0.886)
BLACK = (0.078, 0.078, 0.086)
RED = (0.839, 0.267, 0.184)
BLUE = (0.169, 0.298, 0.608)
YELLOW = (0.914, 0.706, 0.153)
GREY_L = (0.620, 0.606, 0.573)
GREY_D = (0.420, 0.424, 0.447)

DISPLAY, LABEL, MONO = "Century Gothic", "Arial", "Consolas"

SCENES = [  # name, bars, paper background?
    ("title", 4, True), ("problem", 11, False), ("usecase", 12, True),
    ("benefits", 10, False), ("outro", 8, True), ("end", 3, False),
]
STARTS = np.cumsum([0.0] + [b * BAR for _, b, _ in SCENES])
TOTAL = float(STARTS[-1])
TRANSITIONS = ["circle", "columns", "diagonal", "blocks", "bars"]
DATA: dict = {}

# Every schedule the picture and the sound share. Times are seconds inside their own beat.
TIMING = {
    "title": {"circle": 0.25, "square": 0.75, "halfdisc": 1.05, "rule": 1.45, "word": 1.95, "sub": 2.7, "mark": 3.3},
    "p1": {"row0": 0.45, "row_step": 0.045, "row_n": 30, "flip": 2.5, "gene": 3.3, "drop": 4.1},
    "p2": {"c1": 0.4, "c2": 0.8, "arrow": 1.8, "shrink": 2.5, "bad": 4.2},
    "p3": {"seg0": 0.5, "seg_step": 0.22, "seg_n": 12, "grid": 3.8, "grid_step": 0.035, "grid_n": 60, "out": 6.4},
    "u1": {"gene": 0.4, "c1": 1.0, "c2": 1.5, "rule": 2.2},
    "u2": {"bar0": 0.5, "bar_step": 0.035, "bar_n": 44, "pin0": 2.4, "pin_step": 0.18, "pin_n": 6, "count": 4.0},
    "u3": {"row0": 0.5, "row_step": 0.42, "row_n": 5, "note": 3.2},
    "u4": {"sq0": 0.4, "sq_step": 0.3, "sq_n": 6, "verdict": 3.4},
    "b1": {"bar0": 0.5, "bar_step": 0.09, "bar_n": 26, "dot": 3.6, "label": 4.2},
    "b2": {"dots": 0.4, "split": 1.7, "pick": 3.1, "travel": 4.0},
    "b3": {"n0": 0.6, "n_step": 0.055, "n_n": 40, "bar": 3.3, "cand0": 4.1, "cand_step": 0.16},
    "o1": {"type0": 0.5, "type_rate": 26.0, "panel0": 2.6, "panel_step": 0.3, "panel_n": 4},
    "o2": {"type0": 0.4, "type_rate": 30.0, "note0": 2.2, "note_step": 0.45, "note_n": 3},
    "o3": {"card0": 0.4, "card_step": 0.55, "card_n": 3, "item_step": 0.16},
    "end": {"circle": 0.3, "word": 0.9, "rule": 1.5, "pip": 2.1, "url": 3.0},
}


# ---- easing --------------------------------------------------------------------------------------
def clamp(x: float) -> float:
    return max(0.0, min(1.0, x))


def out_cubic(x: float) -> float:
    return 1 - (1 - clamp(x)) ** 3


def in_out(x: float) -> float:
    x = clamp(x)
    return x * x * (3 - 2 * x)


def out_back(x: float, k: float = 1.7) -> float:
    x = clamp(x)
    return 1 + (k + 1) * (x - 1) ** 3 + k * (x - 1) ** 2


def lerp(a: float, b: float, t: float) -> float:
    return a + (b - a) * t


def fade(t: float, t0: float, dur: float = 0.5) -> float:
    return in_out((t - t0) / dur)


# ---- primitives ----------------------------------------------------------------------------------
def text(c, s, x, y, size, color, *, align="left", alpha=1.0, family=LABEL, bold=True,
         dy=0.0, rot=0.0, track=0.0, fit=True):
    """Draw text with its optical centre on y. `track` adds letter spacing, which the small labels want."""
    if alpha <= 0.004 or not s:
        return 0.0
    c.save()
    c.select_font_face(family, cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD if bold else cairo.FONT_WEIGHT_NORMAL)
    c.set_font_size(size)
    widths = [c.text_extents(ch).x_advance + track for ch in s]
    total = sum(widths) - track
    sc = 1.0
    if fit and total > 0:
        avail = {"left": W - x - 80, "right": x - 80}.get(align, 2 * min(x, W - x) - 100)
        if total > avail:
            sc = avail / total
    c.translate(x, y + dy)
    c.rotate(rot)
    c.scale(sc, sc)
    ox = {"left": 0.0, "center": -total / 2, "right": -total}[align]
    c.set_source_rgba(*color, alpha)
    cx = ox
    for ch, adv in zip(s, widths):
        c.move_to(cx, size * 0.35)
        c.show_text(ch)
        cx += adv
    c.restore()
    return total * sc


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


def half_disc(c, x, y, r, color, a0, alpha=1.0):
    """A half disc starting at angle a0: the Bauhaus workhorse."""
    if r <= 0 or alpha <= 0:
        return
    c.new_path()
    c.set_source_rgba(*color, alpha)
    c.move_to(x, y)
    c.arc(x, y, r, a0, a0 + math.pi)
    c.close_path()
    c.fill()


def rect(c, x, y, w, h, color, alpha=1.0, stroke=0.0):
    if w <= 0 or h <= 0 or alpha <= 0:
        return
    c.new_path()
    c.rectangle(x, y, w, h)
    c.set_source_rgba(*color, alpha)
    if stroke:
        c.set_line_width(stroke)
        c.stroke()
    else:
        c.fill()


def rule(c, x0, y0, x1, y1, w, color, alpha=1.0, cap=cairo.LINE_CAP_BUTT):
    if alpha <= 0:
        return
    c.new_path()
    c.set_source_rgba(*color, alpha)
    c.set_line_width(w)
    c.set_line_cap(cap)
    c.move_to(x0, y0)
    c.line_to(x1, y1)
    c.stroke()


def slash(c, x, y, r, color, w, alpha=1.0):
    rule(c, x - r * 0.7, y + r * 0.7, x + r * 0.7, y - r * 0.7, w, color, alpha, cairo.LINE_CAP_ROUND)


# ---- chapter furniture ---------------------------------------------------------------------------
def chapter(c, t, num, label, fg, *, t0=0.2):
    """A small numbered marker that keeps the film's structure visible."""
    a = fade(t, t0, 0.6)
    if a <= 0.004:
        return
    rule(c, W - 420, 112, W - 80, 112, 4, fg, a * 0.9)
    text(c, num, W - 420, 172, 34, fg, alpha=a, family=DISPLAY)
    text(c, label, W - 80, 172, 26, fg, align="right", alpha=a * 0.72, family=LABEL, bold=False, track=5)


def headline(c, t, s, x, y, size, fg, *, t0=0.3, sub=None, sub_gap=76):
    a = fade(t, t0, 0.55)
    dx = (1 - out_cubic((t - t0) / 0.7)) * 36
    text(c, s, x - dx, y, size, fg, alpha=a, family=DISPLAY)
    if sub:
        a2 = fade(t, t0 + 0.4, 0.6)
        text(c, sub, x, y + sub_gap, 32, fg, alpha=a2 * 0.68, family=LABEL, bold=False)


def stage(c, t, beats, fade_dur=0.45, shift=70.0):
    """Run (duration, fn) beats, sliding one out as the next slides in."""
    starts = [0.0]
    for d, _ in beats:
        starts.append(starts[-1] + d)
    for i, (d, fn) in enumerate(beats):
        t_in = starts[i] - (fade_dur if i else 0.0)
        t_out = starts[i + 1]
        if t < t_in or (i < len(beats) - 1 and t >= t_out):
            continue
        a, off = 1.0, 0.0
        if i:
            p = in_out((t - t_in) / fade_dur)
            a, off = p, (1 - p) * shift
        if i < len(beats) - 1:
            p = in_out((t - (t_out - fade_dur)) / fade_dur)
            a, off = min(a, 1 - p), off - p * shift
        if a <= 0.004:
            continue
        c.push_group()
        c.save()
        c.translate(off, 0)
        fn(c, t - t_in)
        c.restore()
        c.pop_group_to_source()
        c.paint_with_alpha(a)


# ---- data ----------------------------------------------------------------------------------------
def load_data() -> dict:
    rep = json.loads((OUT / "report_crc-pd1-exhaustion.json").read_text(encoding="utf-8"))
    its = rep["iterations"]
    null = rep["efficacy_null"]
    cands = []
    for it in its:
        t1 = it["report"]["tier1_selectivity"]
        h = t1.get("healthy_control_deltas") or {}
        hv = next(iter(h.values()), 0.0) if isinstance(h, dict) else 0.0
        cands.append({"pos": it["report"]["candidate"]["target_edit"]["position"],
                      "lfc": t1["target_delta_sketch_b"],
                      "healthy": float(hv) if isinstance(hv, (int, float)) else 0.0})
    meth = [it["report"]["tier2_hazards"]["cpg_methylation_mean_abs_shift"] for it in its]
    over = [v for v in meth if v > 0.12]
    return {"cands": cands, "abs": [abs(x["lfc"]) for x in cands], "null": null["abs_lfc_values"],
            "cutoff": null["abs_lfc_max"], "n_edits": len(its), "n_null": null["n_random_edits"],
            "meth_closest": min(over) if over else max(meth),
            "n_fail_meth": sum(1 for it in its if not it["report"]["tier2_hazards"]["methylation_hazard_passed"])}


# ================================================================================================
# TITLE
# ================================================================================================
def scene_title(c, t, dur, paper):
    T = TIMING["title"]
    p = out_back(fade(t, T["circle"], 0.8), 1.2)
    disc(c, 250, 560, 430 * p, BLUE)
    p = out_cubic(fade(t, T["square"], 0.6))
    rect(c, 1430, lerp(-220, 140, p), 200, 200, YELLOW)
    p = out_cubic(fade(t, T["halfdisc"], 0.7))
    half_disc(c, 1180, H + 40, 250 * p, RED, math.pi)
    p = out_cubic(fade(t, T["rule"], 0.7))
    rule(c, 740, 626, lerp(740, 1800, p), 626, 8, BLACK)
    a = fade(t, T["word"], 0.6)
    dx = (1 - out_cubic((t - T["word"]) / 0.8)) * 44
    w = text(c, "gene-edit-", 740 - dx, 540, 104, BLACK, alpha=a, family=DISPLAY)
    text(c, "check", 740 - dx + w, 540, 104, RED, alpha=a, family=DISPLAY)
    a = fade(t, T["sub"], 0.7)
    text(c, "check a gene edit on a computer,", 744, 700, 36, BLACK, alpha=a * 0.8, family=LABEL, bold=False)
    text(c, "before anyone spends weeks on it in a lab.", 744, 750, 36, BLACK, alpha=a * 0.8, family=LABEL, bold=False)
    a = fade(t, T["mark"], 0.6)
    disc(c, 250, 560, 62, PAPER, a)
    rect(c, 238, 488, 24, 24, RED, a)


# ================================================================================================
# 01  PROBLEM
# ================================================================================================
def _p1(c, t):
    """A gene is a row of letters; change one and the level it sits at changes."""
    T = TIMING["p1"]
    x0, y, s, gap = 150, 500, 30, 14
    for i in range(T["row_n"]):
        p = out_back(fade(t, T["row0"] + i * T["row_step"], 0.28), 2.2)
        if p <= 0:
            continue
        flipped = i == 19 and t > T["flip"]
        sz = s * p * (1.45 if flipped and t < T["flip"] + 0.3 else 1.0)
        rect(c, x0 + i * (s + gap) + (s - sz) / 2, y - sz / 2, sz, sz, RED if flipped else PAPER)
    if t > T["flip"]:
        pr = clamp((t - T["flip"]) / 0.8)
        ring(c, x0 + 19 * (s + gap) + s / 2, y, 30 + 120 * pr, RED, 4, 1 - pr)
    a = fade(t, T["gene"], 0.6)
    r = lerp(170, 74, out_cubic((t - T["drop"]) / 1.4))
    rule(c, 150, 790, 1560 - r - 40, 790, 4, PAPER, a * 0.35)
    disc(c, 1560, 790, r, BLUE, a)
    text(c, "the gene's level", 1560, 990, 28, PAPER, align="center", alpha=a * 0.7, family=LABEL,
         bold=False, track=4)
    headline(c, t, "one letter decides", 150, 260, 96, PAPER, t0=0.15,
             sub="a single base can switch a gene down")


def _p2(c, t):
    """The same edit lands in the cells you wanted to protect."""
    T = TIMING["p2"]
    sh = out_cubic((t - T["shrink"]) / 1.5)
    disc(c, 620, 700, 210 * out_back(fade(t, T["c1"], 0.5)) * lerp(1.0, 0.55, sh), RED)
    disc(c, 1380, 700, 210 * out_back(fade(t, T["c2"], 0.5)) * lerp(1.0, 0.62, sh), BLUE)
    text(c, "cells to change", 620, 980, 28, PAPER, align="center", alpha=fade(t, T["c1"] + 0.3) * 0.7,
         family=LABEL, bold=False, track=4)
    text(c, "cells to protect", 1380, 980, 28, PAPER, align="center", alpha=fade(t, T["c2"] + 0.3) * 0.7,
         family=LABEL, bold=False, track=4)
    p = out_cubic(fade(t, T["arrow"], 0.8))
    if p > 0:
        rule(c, 150, 450, lerp(150, 1700, p), 450, 8, YELLOW)
        rect(c, lerp(150, 1700, p) - 18, 432, 36, 36, YELLOW)
    ab = fade(t, T["bad"], 0.5)
    if ab > 0:
        s = out_back(ab, 2.0) * 150
        ring(c, 1380, 700, s, RED, 10, ab)
        slash(c, 1380, 700, s, RED, 10, ab)
    headline(c, t, "it can hit the wrong cells", 150, 260, 96, PAPER, t0=0.15,
             sub="the ones you needed to leave alone")


def _p3(c, t):
    """Weeks per candidate, times sixty candidates."""
    T = TIMING["p3"]
    y = 520
    for i in range(T["seg_n"]):
        p = out_cubic(fade(t, T["seg0"] + i * T["seg_step"], 0.25))
        if p <= 0:
            continue
        rect(c, 150 + i * 132, y - 34, 112 * p, 68, PAPER if i % 2 else GREY_D)
    weeks = int(clamp((t - T["seg0"]) / (T["seg_step"] * T["seg_n"])) * T["seg_n"])
    text(c, f"{weeks} weeks", 150, y + 128, 46, YELLOW, alpha=fade(t, T["seg0"] + 0.3), family=DISPLAY)
    text(c, "to test one candidate at a bench", 150, y + 186, 28, PAPER,
         alpha=fade(t, T["seg0"] + 0.6) * 0.62, family=LABEL, bold=False, track=3)
    cols, s, gap = 20, 38, 10
    for i in range(T["grid_n"]):
        p = out_back(fade(t, T["grid"] + i * T["grid_step"], 0.22), 2.2)
        if p <= 0:
            continue
        rect(c, 150 + (i % cols) * (s + gap), 790 + (i // cols) * (s + gap), s * p, s * p, RED)
    a = fade(t, T["out"], 0.6)
    text(c, "60", 1230, 860, 150, YELLOW, alpha=a, family=DISPLAY)
    text(c, "candidates", 1230, 950, 30, PAPER, alpha=a * 0.7, family=LABEL, bold=False, track=4)
    headline(c, t, "and the lab takes weeks", 150, 260, 96, PAPER, t0=0.15,
             sub="for every single one of them")


def scene_problem(c, t, dur, paper):
    chapter(c, t, "01", "PROBLEM", PAPER)
    stage(c, t, [(8.0, _p1), (8.2, _p2), (8.2, _p3)])


# ================================================================================================
# 02  USE CASE
# ================================================================================================
def _u1(c, t):
    """What you hand it: one gene, and two cell states."""
    T = TIMING["u1"]
    p = out_back(fade(t, T["gene"], 0.5))
    rect(c, 230, 560 - 110 * p, 220 * p, 220 * p, YELLOW)
    text(c, "a gene", 340, 730, 28, BLACK, align="center", alpha=fade(t, T["gene"] + 0.3) * 0.7,
         family=LABEL, bold=False, track=4)
    disc(c, 1120, 470, 130 * out_back(fade(t, T["c1"], 0.5)), RED)
    disc(c, 1120, 770, 130 * out_back(fade(t, T["c2"], 0.5)), BLUE)
    text(c, "cells to change", 1300, 470, 32, BLACK, alpha=fade(t, T["c1"] + 0.3) * 0.75, family=LABEL, bold=False)
    text(c, "cells to protect", 1300, 770, 32, BLACK, alpha=fade(t, T["c2"] + 0.3) * 0.75, family=LABEL, bold=False)
    pr = out_cubic(fade(t, T["rule"], 0.8))
    if pr > 0:
        rule(c, 470, 560, lerp(470, 930, pr), 560, 6, BLACK)
        rule(c, 930, 560, 930, lerp(560, 470, pr), 6, BLACK)
        rule(c, 930, 560, 930, lerp(560, 770, pr), 6, BLACK)
        rule(c, 930, 470, lerp(930, 985, pr), 470, 6, BLACK)
        rule(c, 930, 770, lerp(930, 985, pr), 770, 6, BLACK)
    headline(c, t, "you give it two things", 150, 260, 92, BLACK, t0=0.15,
             sub="a gene, and the two kinds of cell you care about")


_U2_H = [0.14, 0.2, 0.33, 0.25, 0.18, 0.4, 0.62, 0.48, 0.3, 0.22, 0.35, 0.55, 0.78, 0.95, 0.72, 0.5,
         0.36, 0.28, 0.42, 0.6, 0.47, 0.33, 0.25, 0.4, 0.58, 0.44, 0.3, 0.2, 0.28, 0.46, 0.68, 0.52,
         0.35, 0.24, 0.3, 0.44, 0.33, 0.22, 0.28, 0.38, 0.26, 0.18, 0.22, 0.15]
_U2_PINS = [6, 12, 13, 19, 30, 36]


def _u2(c, t):
    """Where the model is sensitive, and what an editor can actually reach."""
    T = TIMING["u2"]
    base, x0, bw = 840, 170, 36
    for i in range(T["bar_n"]):
        p = out_cubic(fade(t, T["bar0"] + i * T["bar_step"], 0.3))
        if p <= 0:
            continue
        h = _U2_H[i] * 420 * p
        rect(c, x0 + i * (bw + 4), base - h, bw, h, BLUE if _U2_H[i] >= 0.6 else BLACK)
    for k in range(T["pin_n"]):
        p = out_back(fade(t, T["pin0"] + k * T["pin_step"], 0.3), 2.2)
        if p <= 0:
            continue
        i = _U2_PINS[k]
        x = x0 + i * (bw + 4) + bw / 2
        h = _U2_H[i] * 420
        rule(c, x, base - h - 10, x, base - h - 58, 4, RED, p)
        disc(c, x, base - h - 76, 15 * p, RED)
    a = fade(t, T["count"], 0.6)
    text(c, f"{int(60 * out_cubic((t - T['count']) / 1.0))}", 1810, 350, 130, RED, align="right",
         alpha=a, family=DISPLAY)
    text(c, "edits it could make", 1810, 436, 28, BLACK, align="right", alpha=a * 0.7,
         family=LABEL, bold=False, track=4)
    headline(c, t, "it finds every edit", 150, 260, 92, BLACK, t0=0.15,
             sub="ranked by where the model is most sensitive")


def _u3(c, t):
    """Each candidate becomes two marks: how far the target moved, how far the others did."""
    T = TIMING["u3"]
    mid, scale = 1520, 3300.0
    rule(c, mid, 430, mid, 900, 4, BLACK, fade(t, 0.2) * 0.35)
    text(c, "no change", mid, 946, 24, BLACK, align="center", alpha=fade(t, 0.3) * 0.5,
         family=LABEL, bold=False, track=4)
    for i, cd in enumerate(DATA["cands"][:T["row_n"]]):
        y = 470 + i * 96
        p = out_cubic(fade(t, T["row0"] + i * T["row_step"], 0.55))
        if p <= 0:
            continue
        wt = abs(2 ** cd["lfc"] - 1) * scale * p
        wh = abs(2 ** cd["healthy"] - 1) * scale * p
        rect(c, mid - wt, y - 32, wt, 28, RED)
        rect(c, mid - wh, y + 2, wh, 28, BLUE)
        if wt - wh > 60:
            rect(c, mid + 34, y - 11, 22, 22, BLACK, fade(t, T["row0"] + i * T["row_step"] + 0.5))
    a = fade(t, T["note"], 0.6)
    rect(c, 200, 492, 34, 28, RED, a)
    text(c, "cells to change", 252, 506, 28, BLACK, alpha=a * 0.8, family=LABEL, bold=False)
    rect(c, 200, 548, 34, 28, BLUE, a)
    text(c, "cells to protect", 252, 562, 28, BLACK, alpha=a * 0.8, family=LABEL, bold=False)
    text(c, "a wide gap is the whole point", 200, 650, 28, BLACK,
         alpha=fade(t, T["note"] + 0.6) * 0.6, family=LABEL, bold=False)
    headline(c, t, "it predicts each one", 150, 260, 92, BLACK, t0=0.15,
             sub="in both kinds of cell at once")


def _u4(c, t):
    """Six limits, then a verdict."""
    T = TIMING["u4"]
    s, gap = 190, 44
    x0 = (W - (T["sq_n"] * s + (T["sq_n"] - 1) * gap)) / 2
    for i in range(T["sq_n"]):
        p = out_back(fade(t, T["sq0"] + i * T["sq_step"], 0.3), 2.0)
        if p <= 0:
            continue
        bad = i == 3
        x, y = x0 + i * (s + gap), 440
        rect(c, x, y, s, s, BLACK, p * 0.1)
        fill = out_cubic(fade(t, T["sq0"] + i * T["sq_step"] + 0.25, 0.4))
        if fill > 0:
            rect(c, x, y + s * (1 - fill), s, s * fill, RED if bad else BLUE)
        if bad and fill > 0.9:
            slash(c, x + s / 2, y + s / 2, s * 0.3, PAPER, 10, (fill - 0.9) * 10)
    a = fade(t, T["verdict"], 0.7)
    if a > 0:
        p = out_back(a, 1.6) * 140
        disc(c, 1600, 860, p, RED, a)
        slash(c, 1600, 860, p, PAPER, 12, a)
        text(c, "none passed", 1420, 850, 46, BLACK, align="right", alpha=a, family=DISPLAY)
        text(c, f"{DATA['n_fail_meth']} of {DATA['n_edits']} stopped by one check", 1420, 906, 28, BLACK,
             align="right", alpha=a * 0.7, family=LABEL, bold=False)
    headline(c, t, "then it tries to break them", 150, 260, 92, BLACK, t0=0.15,
             sub="six limits, every one fixed before the run")


def scene_usecase(c, t, dur, paper):
    chapter(c, t, "02", "USE CASE", BLACK)
    stage(c, t, [(6.6, _u1), (6.7, _u2), (6.7, _u3), (6.6, _u4)])


# ================================================================================================
# 03  BENEFITS
# ================================================================================================
def _b1(c, t):
    """Bench time against one run."""
    T = TIMING["b1"]
    y = 540
    for i in range(T["bar_n"]):
        p = out_cubic(fade(t, T["bar0"] + i * T["bar_step"], 0.22))
        if p <= 0:
            continue
        rect(c, 150 + i * 64, y - 40, 52 * p, 80 * p, PAPER if i % 2 else GREY_D)
    text(c, "the lab", 150, y + 110, 30, PAPER, alpha=fade(t, T["bar0"] + 0.4) * 0.7,
         family=LABEL, bold=False, track=4)
    a = out_back(fade(t, T["dot"], 0.5), 2.2)
    if a > 0:
        rect(c, 150, y + 250, 52 * a, 80 * a, YELLOW)
        text(c, "one run", 150, y + 400, 30, YELLOW, alpha=fade(t, T["dot"] + 0.2) * 0.9,
             family=LABEL, bold=False, track=4)
    text(c, "minutes", 1800, y + 300, 112, YELLOW, align="right", alpha=fade(t, T["label"], 0.7), family=DISPLAY)
    headline(c, t, "weeks become minutes", 150, 260, 96, PAPER, t0=0.15,
             sub="the obvious losers never reach a bench")


def _b2(c, t):
    """Propose on one half of the patients, judge on the other."""
    T = TIMING["b2"]
    rng = np.random.default_rng(21)
    pts = [(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(44)]
    sp = out_cubic((t - T["split"]) / 1.6)
    for i, (ux, uy) in enumerate(pts):
        a = fade(t, T["dots"] + i * 0.016, 0.4)
        grp = i % 2
        cx = lerp(960, 580 if grp == 0 else 1380, sp)
        col = PAPER if sp < 0.3 else (RED if grp == 0 else BLUE)
        rect(c, cx + ux * lerp(330, 190, sp) - 11, 650 + uy * lerp(240, 190, sp) - 11, 22, 22, col, a)
    if sp > 0.05:
        rule(c, 960, 370, 960, 950, 6, YELLOW, sp)
        text(c, "split by patient", 960, 334, 26, YELLOW, align="center", alpha=sp,
             family=LABEL, bold=False, track=4)
    pick = out_back(fade(t, T["pick"], 0.4), 2.0)
    if pick > 0:
        rect(c, 580 - 26 * pick, 650 - 26 * pick, 52 * pick, 52 * pick, YELLOW)
    tr = out_cubic((t - T["travel"]) / 1.4)
    if tr > 0:
        x = lerp(580, 1380, tr)
        rule(c, 580, 650, x, 650, 4, YELLOW, 0.6)
        rect(c, x - 26, 624, 52, 52, YELLOW)
        if tr > 0.95:
            ring(c, 1380, 650, 190, YELLOW, 5, clamp((tr - 0.95) * 20))
    headline(c, t, "judged on patients it never saw", 150, 260, 86, PAPER, t0=0.15,
             sub="one half proposes, the other half decides")


def _b3(c, t):
    """Nothing counts until it clears a crowd of random edits."""
    T = TIMING["b3"]
    ax0, ax1, ay = 220, 1700, 820
    rule(c, ax0, ay, ax1, ay, 4, PAPER, fade(t, 0.3) * 0.5)
    xmax = max(DATA["abs"]) * 1.12
    jit = np.random.default_rng(77).normal(0, 13, len(DATA["null"]))

    def px(v):
        return ax0 + 40 + math.sqrt(min(v, xmax) / xmax) * (ax1 - ax0 - 80)

    for i, v in enumerate(DATA["null"]):
        p = out_back(fade(t, T["n0"] + i * T["n_step"], 0.25), 2.2)
        if p <= 0:
            continue
        rect(c, px(v) + jit[i] - 11 * p, ay - 40 - 34 * (i % 5) - 11 * p, 22 * p, 22 * p, GREY_L)
    text(c, f"{DATA['n_null']} random edits", ax0, 600, 32, GREY_L, alpha=fade(t, T["n0"] + 0.5) * 0.9,
         family=LABEL, bold=False, track=3)
    a = fade(t, T["bar"], 0.6)
    cx = px(DATA["cutoff"])
    rule(c, cx, ay + 40, cx, ay - 300, 6, YELLOW, a)
    text(c, "the bar", cx + 24, ay - 326, 30, YELLOW, alpha=a, family=LABEL, bold=False, track=3)
    for i, v in enumerate(DATA["abs"]):
        p = out_back(fade(t, T["cand0"] + i * T["cand_step"], 0.3), 2.0)
        if p <= 0:
            continue
        s = 34 * p
        rect(c, px(v) - s / 2, ay - 60 - 50 * (i % 3) - s / 2, s, s, RED if v > DATA["cutoff"] else GREY_D)
    headline(c, t, "and it has to beat chance", 150, 260, 96, PAPER, t0=0.15,
             sub="a crowd of random edits sets the bar first")


def scene_benefits(c, t, dur, paper):
    chapter(c, t, "03", "WHAT YOU GET", PAPER)
    stage(c, t, [(7.4, _b1), (7.4, _b2), (7.4, _b3)])


# ================================================================================================
# 04  OUTRO
# ================================================================================================
def _o1(c, t):
    T = TIMING["o1"]
    half_disc(c, 1700, H + 30, 300 * out_cubic(fade(t, 0.3, 0.9)), BLUE, math.pi)
    rect(c, 150, 430, 940, 290, BLACK, fade(t, 0.15) * 0.06)
    rule(c, 150, 430, 150, 720, 8, BLACK, fade(t, 0.15))
    shown = int(max(0.0, t - T["type0"]) * T["type_rate"])
    y = 510
    for s in ("$ pip install gene-edit-check", "$ gene-edit-check doctor"):
        n = min(len(s), shown)
        shown -= max(n, 0)
        if n <= 0:
            break
        text(c, s[:n], 200, y, 38, BLACK, family=MONO, bold=False, fit=False)
        y += 80
    for i, (k, col) in enumerate((("api key", RED), ("reference data", BLUE), ("cache", BLUE), ("ready", YELLOW))):
        p = out_back(fade(t, T["panel0"] + i * T["panel_step"], 0.3), 2.0)
        if p <= 0:
            continue
        y2 = 440 + i * 86
        rect(c, 1250, y2, 54 * p, 54 * p, col)
        text(c, k, 1330, y2 + 27, 32, BLACK, alpha=clamp(p) * 0.85, family=LABEL, bold=False)
    headline(c, t, "install it", 150, 260, 92, BLACK, t0=0.1, sub="one command, then it checks its own setup")


def _o2(c, t):
    T = TIMING["o2"]
    disc(c, 240, H + 120, 300 * out_cubic(fade(t, 0.3, 0.9)), RED)
    rect(c, 150, 420, 1060, 360, BLACK, fade(t, 0.15) * 0.06)
    rule(c, 150, 420, 150, 780, 8, BLACK, fade(t, 0.15))
    shown = int(max(0.0, t - T["type0"]) * T["type_rate"])
    y = 490
    for s in ("$ gene-edit-check run-h5ad \\", "   --h5ad my_cells.h5ad \\",
              "   --gene-symbol PDCD1 \\", "   --target-query \"state == 'exhausted'\""):
        n = min(len(s), shown)
        shown -= max(n, 0)
        if n <= 0:
            break
        text(c, s[:n], 200, y, 34, BLACK, family=MONO, bold=False, fit=False)
        y += 76
    for i, (s, col) in enumerate((("your cells", YELLOW), ("any gene", RED), ("which cells to change", BLUE))):
        p = out_cubic(fade(t, T["note0"] + i * T["note_step"], 0.5))
        if p <= 0:
            continue
        y2 = 566 + i * 76
        rule(c, 1250, y2, lerp(1250, 1330, p), y2, 5, col)
        text(c, s, 1350, y2, 30, BLACK, alpha=p * 0.85, family=LABEL, bold=False)
    headline(c, t, "point it at your cells", 150, 260, 92, BLACK, t0=0.1,
             sub="any gene, any single-cell dataset")


def _o3(c, t):
    T = TIMING["o3"]
    cards = [("a report", ["every candidate", "every check", "the thresholds"], RED),
             ("the figures", None, BLUE),
             ("a candidate", ["position", "editor and guide", "atlas link"], YELLOW)]
    cw, gap = 490, 60
    x0 = (W - (3 * cw + 2 * gap)) / 2
    for i, (title_, rows, col) in enumerate(cards):
        p = out_cubic(fade(t, T["card0"] + i * T["card_step"], 0.5))
        if p <= 0:
            continue
        x, y = x0 + i * (cw + gap), 440 + (1 - p) * 60
        rect(c, x, y, cw, 70, col, p)
        text(c, title_, x + 28, y + 35, 34, PAPER if i != 2 else BLACK, alpha=p, family=DISPLAY)
        rect(c, x, y + 70, cw, 310, BLACK, p * 0.05)
        if rows:
            for j, r in enumerate(rows):
                a = fade(t, T["card0"] + i * T["card_step"] + 0.3 + j * T["item_step"], 0.35)
                rect(c, x + 28, y + 132 + j * 74 - 9, 18, 18, col, a)
                text(c, r, x + 64, y + 132 + j * 74, 30, BLACK, alpha=a * 0.85, family=LABEL, bold=False)
        else:
            for j in range(3):
                for k in range(4):
                    a = fade(t, T["card0"] + i * T["card_step"] + 0.3 + (j * 4 + k) * 0.05, 0.3)
                    rect(c, x + 46 + k * 102, y + 124 + j * 88, 68, 68,
                         RED if (k == 2 and j != 1) else BLUE, a * 0.9)
    headline(c, t, "it writes the evidence down", 150, 260, 92, BLACK, t0=0.1,
             sub="not a verdict you have to take on trust")


def scene_outro(c, t, dur, paper):
    chapter(c, t, "04", "GET IT", BLACK)
    stage(c, t, [(5.9, _o1), (5.9, _o2), (6.0, _o3)])


def scene_end(c, t, dur, paper):
    T = TIMING["end"]
    disc(c, 1640, 300, 300 * out_back(fade(t, T["circle"], 0.7), 1.3), BLUE)
    half_disc(c, 1320, H, 210 * out_cubic(fade(t, T["circle"] + 0.3, 0.7)), YELLOW, math.pi)
    a = fade(t, T["word"], 0.6)
    dx = (1 - out_cubic((t - T["word"]) / 0.8)) * 40
    w = text(c, "gene-edit-", 150 - dx, 470, 108, PAPER, alpha=a, family=DISPLAY)
    text(c, "check", 150 - dx + w, 470, 108, RED, alpha=a, family=DISPLAY)
    rule(c, 150, 560, lerp(150, 1180, out_cubic(fade(t, T["rule"], 0.7))), 560, 8, PAPER)
    text(c, "pip install gene-edit-check", 150, 652, 44, YELLOW, alpha=fade(t, T["pip"], 0.6),
         family=MONO, bold=False)
    a = fade(t, T["url"], 0.6)
    text(c, "github.com/aravmendiratta/gene-edit-check", 150, 742, 32, PAPER, alpha=a * 0.75,
         family=LABEL, bold=False)
    for i, col in enumerate((RED, BLUE, YELLOW)):
        rect(c, 150 + i * 56, 840, 40, 40, col, a)


SCENE_FN = {"title": scene_title, "problem": scene_problem, "usecase": scene_usecase,
            "benefits": scene_benefits, "outro": scene_outro, "end": scene_end}


# ---- transitions ---------------------------------------------------------------------------------
def clip_shape(c, kind, p):
    """Clip to the shape that reveals the incoming chapter."""
    c.new_path()
    if kind == "circle":
        c.arc(1420, 360, math.hypot(W, H) * 1.05 * in_out(p), 0, 2 * math.pi)
    elif kind == "columns":
        n = 9
        for i in range(n):
            h = H * in_out(clamp((p - (i % 3) * 0.09) / 0.72))
            c.rectangle(i * W / n, 0 if i % 2 == 0 else H - h, W / n + 1, h)
    elif kind == "diagonal":
        q = in_out(p) * (W + H) * 1.1
        c.move_to(-H, 0)
        c.line_to(-H + q, 0)
        c.line_to(q, H)
        c.line_to(-H, H)
        c.close_path()
    elif kind == "blocks":
        cols, rows_ = 8, 5
        bw, bh = W / cols, H / rows_
        for i in range(cols):
            for j in range(rows_):
                s = out_back(clamp((p - ((i + j) % 5) * 0.07) / 0.68), 1.4)
                if s <= 0:
                    continue
                c.rectangle(i * bw + bw * (1 - s) / 2, j * bh + bh * (1 - s) / 2, bw * s + 1, bh * s + 1)
    else:  # bars
        n = 7
        for i in range(n):
            w = W * in_out(clamp((p - i * 0.06) / 0.7))
            c.rectangle(0 if i % 2 == 0 else W - w, i * H / n, w, H / n + 1)
    c.clip()


def draw_scene(c, idx, t):
    name, bars, paper = SCENES[idx]
    c.set_source_rgb(*(PAPER if paper else BLACK))
    c.paint()
    SCENE_FN[name](c, t, bars * BAR, paper)


TRANS_DUR = 0.75


def render(t, surf):
    c = cairo.Context(surf)
    c.scale(OW / W, OH / H)
    idx = max(i for i in range(len(SCENES)) if STARTS[i] <= t + 1e-9)
    local = t - STARTS[idx]
    if idx > 0 and local < TRANS_DUR:
        draw_scene(c, idx - 1, SCENES[idx - 1][1] * BAR + local)
        c.save()
        clip_shape(c, TRANSITIONS[(idx - 1) % len(TRANSITIONS)], local / TRANS_DUR)
        draw_scene(c, idx, local)
        c.restore()
    else:
        draw_scene(c, idx, local)


# ================================================================================================
# sound
# ================================================================================================
SR = 44100
SCALE = [220.00, 261.63, 293.66, 329.63, 392.00, 440.00, 523.25, 587.33, 659.25]  # A minor


def build_events() -> list[tuple[float, str, float]]:
    """Every sound effect, mirroring the animation schedules in TIMING."""
    ev: list[tuple[float, str, float]] = []
    S = {name: float(STARTS[i]) for i, (name, _, _) in enumerate(SCENES)}

    def add(scene, t, kind, pitch=0.0):
        ev.append((S[scene] + t, kind, pitch))

    T = TIMING["title"]
    for k, key in enumerate(("circle", "square", "halfdisc")):
        add("title", T[key], "thunk", SCALE[k])
    add("title", T["rule"], "whoosh")
    add("title", T["word"], "chime", SCALE[5])
    add("title", T["sub"], "bop", SCALE[2])
    add("title", T["mark"], "beep", SCALE[6])

    T = TIMING["p1"]
    for i in range(T["row_n"]):
        add("problem", T["row0"] + i * T["row_step"], "tick", SCALE[i % 4] * 2)
    add("problem", T["flip"], "beep", SCALE[6])
    add("problem", T["gene"], "bop", SCALE[3])
    add("problem", T["drop"], "fall", SCALE[1])
    b, T = 8.0, TIMING["p2"]
    add("problem", b + T["c1"], "thunk", SCALE[2])
    add("problem", b + T["c2"], "thunk", SCALE[0])
    add("problem", b + T["arrow"], "whoosh")
    add("problem", b + T["shrink"], "fall", SCALE[2])
    add("problem", b + T["bad"], "alarm", SCALE[1])
    b, T = 16.2, TIMING["p3"]
    for i in range(T["seg_n"]):
        add("problem", b + T["seg0"] + i * T["seg_step"], "bop", SCALE[i % 5] / 2)
    for i in range(0, T["grid_n"], 4):
        add("problem", b + T["grid"] + i * T["grid_step"], "tick", SCALE[(i // 4) % 6] * 2)
    add("problem", b + T["out"], "chime", SCALE[4])

    T = TIMING["u1"]
    add("usecase", T["gene"], "thunk", SCALE[4])
    add("usecase", T["c1"], "thunk", SCALE[2])
    add("usecase", T["c2"], "thunk", SCALE[0])
    add("usecase", T["rule"], "whoosh")
    b, T = 6.6, TIMING["u2"]
    for i in range(0, T["bar_n"], 2):
        add("usecase", b + T["bar0"] + i * T["bar_step"], "tick", SCALE[i % 7] * 2)
    for k in range(T["pin_n"]):
        add("usecase", b + T["pin0"] + k * T["pin_step"], "beep", SCALE[k % 6] * 2)
    add("usecase", b + T["count"], "chime", SCALE[5])
    b, T = 13.3, TIMING["u3"]
    for i in range(T["row_n"]):
        add("usecase", b + T["row0"] + i * T["row_step"], "swipe", SCALE[i % 5])
    add("usecase", b + T["note"], "bop", SCALE[3])
    b, T = 20.0, TIMING["u4"]
    for i in range(T["sq_n"]):
        add("usecase", b + T["sq0"] + i * T["sq_step"], "bop", SCALE[i % 6])
        add("usecase", b + T["sq0"] + i * T["sq_step"] + 0.25, "fill", SCALE[i % 6] * 2)
    add("usecase", b + T["verdict"], "alarm", SCALE[0])

    T = TIMING["b1"]
    for i in range(0, T["bar_n"], 2):
        add("benefits", T["bar0"] + i * T["bar_step"], "tick", SCALE[i % 6] * 2)
    add("benefits", T["dot"], "thunk", SCALE[4])
    add("benefits", T["label"], "chime", SCALE[6])
    b, T = 7.4, TIMING["b2"]
    for i in range(0, 44, 4):
        add("benefits", b + T["dots"] + i * 0.016, "tick", SCALE[(i // 4) % 7] * 2)
    add("benefits", b + T["split"], "whoosh")
    add("benefits", b + T["pick"], "beep", SCALE[5])
    add("benefits", b + T["travel"], "swipe", SCALE[3])
    add("benefits", b + T["travel"] + 1.33, "chime", SCALE[6])
    b, T = 14.8, TIMING["b3"]
    for i in range(0, T["n_n"], 3):
        add("benefits", b + T["n0"] + i * T["n_step"], "tick", SCALE[(i // 3) % 7] * 2)
    add("benefits", b + T["bar"], "thunk", SCALE[4])
    for i in range(8):
        add("benefits", b + T["cand0"] + i * T["cand_step"], "bop", SCALE[i % 6])

    T = TIMING["o1"]
    for i in range(14):
        add("outro", T["type0"] + i * 0.12, "key", SCALE[i % 5] * 2)
    for i in range(T["panel_n"]):
        add("outro", T["panel0"] + i * T["panel_step"], "bop", SCALE[i % 5])
    b, T = 5.9, TIMING["o2"]
    for i in range(16):
        add("outro", b + T["type0"] + i * 0.11, "key", SCALE[i % 5] * 2)
    for i in range(T["note_n"]):
        add("outro", b + T["note0"] + i * T["note_step"], "swipe", SCALE[i % 5])
    b, T = 11.8, TIMING["o3"]
    for i in range(T["card_n"]):
        add("outro", b + T["card0"] + i * T["card_step"], "thunk", SCALE[i % 5])
        for j in range(3):
            add("outro", b + T["card0"] + i * T["card_step"] + 0.3 + j * T["item_step"], "tick", SCALE[(i + j) % 7] * 2)

    T = TIMING["end"]
    add("end", T["circle"], "thunk", SCALE[0])
    add("end", T["word"], "chime", SCALE[5])
    add("end", T["rule"], "whoosh")
    add("end", T["pip"], "key", SCALE[6])
    add("end", T["url"], "bop", SCALE[4])

    # the transition sweep, and the hit the new chapter lands on
    for i in range(1, len(SCENES)):
        t0 = float(STARTS[i])
        ev.append((t0 - TRANS_DUR, "sweep", 0.0))
        ev.append((t0, "impact", 0.0))
    return [e for e in ev if 0 <= e[0] < TOTAL]


def make_music() -> np.ndarray:
    """A 108 BPM bed (bass, woodblock pulse, marimba figure) with the sound-effect layer on top."""
    from scipy.signal import butter, fftconvolve, lfilter

    n = int(TOTAL * SR) + SR
    dry = np.zeros((n, 2))
    wet = np.zeros((n, 2))
    rng = np.random.default_rng(5)

    def tt(d):
        return np.arange(max(int(d * SR), 1)) / SR

    def filt(x, kind, fc, order=2):
        b, a = butter(order, min(fc / (SR / 2), 0.99), btype=kind)
        return lfilter(b, a, x)

    def put(buf, t0, sig, gl=1.0, gr=None):
        i = int(t0 * SR)
        if i < 0 or i >= n:
            return
        m = min(len(sig), n - i)
        buf[i : i + m, 0] += sig[:m] * gl
        buf[i : i + m, 1] += sig[:m] * (gl if gr is None else gr)

    def v_tick(f):
        t = tt(0.045)
        return (np.sin(2 * np.pi * max(f, 200) * t) * 0.5 + rng.normal(0, 1, len(t)) * 0.3) * np.exp(-t * 110) * 0.15

    def v_bop(f):
        t = tt(0.16)
        return np.sin(2 * np.pi * f * t) * np.exp(-t * 26) * 0.24

    def v_beep(f):
        t = tt(0.11)
        s = np.sign(np.sin(2 * np.pi * f * t)) * 0.35 + np.sin(2 * np.pi * f * 2 * t) * 0.3
        return filt(s, "low", 4200) * np.exp(-t * 34) * 0.18

    def v_key(f):
        t = tt(0.06)
        return (np.sin(2 * np.pi * f * t) * 0.4 + rng.normal(0, 1, len(t)) * 0.5) * np.exp(-t * 90) * 0.11

    def v_thunk(f):
        t = tt(0.4)
        body = np.sin(2 * np.pi * np.cumsum(max(f, 60) * 0.5 + 90 * np.exp(-t * 26)) / SR) * np.exp(-t * 11)
        return (body + rng.normal(0, 1, len(t)) * np.exp(-t * 150) * 0.4) * 0.3

    def v_fill(f):
        t = tt(0.3)
        return np.sin(2 * np.pi * np.cumsum(f * 0.6 + f * 0.5 * (t / 0.3)) / SR) * np.exp(-t * 9) * 0.1

    def v_fall(f):
        t = tt(0.55)
        return np.sin(2 * np.pi * np.cumsum(f * np.exp(-t * 2.4) + 70) / SR) * np.exp(-t * 4.5) * 0.18

    def v_swipe(f):
        t = tt(0.3)
        return (filt(rng.normal(0, 1, len(t)), "high", 2200) * np.exp(-t * 13) * 0.1
                + np.sin(2 * np.pi * f * t) * np.exp(-t * 16) * 0.05)

    def v_whoosh(_f):
        d = 0.6
        t = tt(d)
        return filt(rng.normal(0, 1, len(t)), "high", 900) * (np.sin(np.pi * t / d) ** 2) * 0.11

    def v_chime(f):
        t = tt(1.6)
        s = (np.sin(2 * np.pi * f * t) + 0.4 * np.sin(2 * np.pi * 2.01 * f * t) * np.exp(-t * 4)
             + 0.2 * np.sin(2 * np.pi * 3 * f * t) * np.exp(-t * 7))
        return s * np.exp(-t * 2.6) * 0.2

    def v_alarm(f):
        t = tt(0.5)
        s = np.sign(np.sin(2 * np.pi * f * t * (1 + 0.04 * np.sin(2 * np.pi * 9 * t))))
        return filt(s, "low", 2000) * np.exp(-t * 7) * 0.15

    def v_sweep(_f):
        d = TRANS_DUR
        t = tt(d)
        return (np.sin(2 * np.pi * np.cumsum(150 + 2400 * (t / d) ** 2) / SR) * (t / d) ** 2 * 0.15
                + filt(rng.normal(0, 1, len(t)), "high", 1400) * (t / d) ** 3 * 0.09)

    def v_impact(_f):
        t = tt(1.2)
        return (np.sin(2 * np.pi * np.cumsum(60 * np.exp(-t * 3) + 34) / SR) * np.exp(-t * 3.4) * 0.4
                + filt(rng.normal(0, 1, len(t)), "high", 2600) * np.exp(-t * 9) * 0.15)

    VOICES = {"tick": v_tick, "bop": v_bop, "beep": v_beep, "key": v_key, "thunk": v_thunk, "fill": v_fill,
              "fall": v_fall, "swipe": v_swipe, "whoosh": v_whoosh, "chime": v_chime, "alarm": v_alarm,
              "sweep": v_sweep, "impact": v_impact}

    chords = [(110.0, [220.0, 261.63, 329.63]), (82.41, [164.81, 220.0, 261.63]),
              (98.0, [196.0, 246.94, 293.66]), (87.31, [174.61, 220.0, 261.63])]

    def woodblock(f):
        t = tt(0.09)
        return (np.sin(2 * np.pi * f * t) * 0.6 + rng.normal(0, 1, len(t)) * 0.25) * np.exp(-t * 60) * 0.08

    def marimba(f, d=0.45):
        t = tt(d)
        s = np.sin(2 * np.pi * f * t) + 0.3 * np.sin(2 * np.pi * 4 * f * t) * np.exp(-t * 18)
        return s * np.exp(-t * 7.5) * 0.075

    for bar in range(int(round(TOTAL / BAR)) + 1):
        root, triad = chords[bar % 4]
        t_bar = bar * BAR
        d = BAR + 0.9
        env = np.minimum(1, tt(d) / 0.7) * np.minimum(1, (d - tt(d)) / 0.9)
        pad_l = sum(2 * ((tt(d) * f * (1 + de)) % 1.0) - 1 for f in triad for de in (-0.003, 0.002)) / 6
        pad_r = sum(2 * ((tt(d) * f * (1 + de)) % 1.0) - 1 for f in triad for de in (0.003, -0.002)) / 6
        put(wet, t_bar, filt(pad_l, "low", 1000) * env * 0.4, 1.0, 0.0)
        put(wet, t_bar, filt(pad_r, "low", 1000) * env * 0.4, 0.0, 1.0)
        put(dry, t_bar, np.sin(2 * np.pi * root * tt(BAR * 1.02)) * np.exp(-tt(BAR * 1.02) * 1.5) * 0.2)
        if bar >= 1:
            for b4 in range(4):
                put(dry, t_bar + b4 * BEAT, woodblock(1400 if b4 == 0 else 1000), 0.7, 0.85)
                put(dry, t_bar + b4 * BEAT + BEAT * 0.5, woodblock(820), 0.5, 0.38)
        if bar >= 2:
            fig = [0, 2, 1, 2, 0, 1, 2, 1]
            for e in range(8):
                f = triad[fig[e] % 3] * (2 if e in (2, 6) else 1)
                pan = 0.4 + 0.22 * (e % 2)
                put(dry, t_bar + e * BEAT / 2, marimba(f), 1 - pan, pan)

    for t0, kind, pitch in build_events():
        sig = VOICES[kind](pitch if pitch else 220.0)
        if kind in ("impact", "sweep", "thunk"):
            put(dry, t0, sig)
        else:
            pan = 0.5 + 0.22 * math.sin(t0 * 2.7)
            put(dry, t0, sig, 1 - pan, pan)

    ir_t = tt(2.2)
    ir = filt(rng.normal(0, 1, len(ir_t)), "low", 3400) * np.exp(-ir_t * 3.0)
    ir /= np.abs(ir).sum() / 11
    rev = np.stack([fftconvolve(wet[:, 0], ir)[:n], fftconvolve(wet[:, 1], np.roll(ir, 173))[:n]], axis=1)
    mix = (dry + rev * 0.85 + wet * 0.25)[: int(TOTAL * SR)]
    m = len(mix)
    mix *= (np.minimum(1, np.arange(m) / SR / 0.5) * np.minimum(1, (m - np.arange(m)) / SR / 1.8))[:, None]
    mix = np.tanh(mix * 0.95)
    mix /= max(1e-6, np.abs(mix).max()) / 0.82
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
           "-c:v", "libx264", "-profile:v", "main", "-level", "4.0", "-preset", "slow", "-crf", "21",
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
