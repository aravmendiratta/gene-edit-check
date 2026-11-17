"""The gene-edit-check launch film: a Neo-Bauhaus explainer, drawn and scored entirely in code.

    pip install pycairo imageio-ffmpeg numpy scipy
    PYTHONPATH=src python scripts/make_launch_video.py            # writes assets/launch.mp4 (1280x720)
    PYTHONPATH=src python scripts/make_launch_video.py --preview out.png 12 34 56

Design: a Bauhaus vocabulary (circle, square, rule, primary red/blue/yellow on paper or near-black) used to
draw real biology: a turning double helix, strands that unwind, cells with membranes and nuclei. Compositions
are asymmetric, shapes bleed off the frame, and each chapter arrives through a different geometric transition.

Structure: title, context (why a single-letter edit is worth checking at all), problem, how it works, what you
get, how to get it. A faint scrolling strand runs under the middle chapters as connective tissue, arrows and
leader lines say what every mark means, and headlines name their subject instead of leaning on "it".

Sound: a 108 BPM bed plus a sound effect on every element that lands, all synthesised here. `build_events()`
and the draw calls both read the schedules in `TIMING` and the beat table in `BEATS`, so the picture and the
sound cannot drift apart.
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
    ("title", 4, True), ("context", 8, False), ("problem", 7, True), ("work", 10, False),
    ("how", 8, False), ("gain", 8, True), ("outro", 4, False), ("end", 3, True),
]
STARTS = np.cumsum([0.0] + [b * BAR for _, b, _ in SCENES])
TOTAL = float(STARTS[-1])
TRANSITIONS = ["circle", "columns", "diagonal", "blocks", "bars", "circle", "columns"]
SPINE_SCENES = {"context", "problem", "work", "how"}  # the chapters the connective strand runs under
DATA: dict = {}

# Every schedule the picture and the sound share. Times are seconds inside their own beat.
TIMING = {
    "title": {"helix": 0.15, "square": 0.7, "halfdisc": 0.95, "rule": 1.25, "word": 1.6, "sub": 2.3, "cells": 2.9},
    "c1": {"helix": 0.3, "bud0": 1.5, "bud_step": 0.35, "bud_n": 4, "read": 3.3},
    "c2": {"strand": 0.3, "editor": 1.4, "flip": 2.4, "label": 3.2},
    "c3": {"cell1": 0.4, "cell2": 0.9, "dim": 1.9, "bad": 3.3},
    "p1": {"strand": 0.3, "pin0": 1.1, "pin_step": 0.16, "pin_n": 8, "count": 2.8},
    "p2": {"seg0": 0.35, "seg_step": 0.14, "seg_n": 12, "grid": 2.4, "grid_step": 0.022, "grid_n": 60},
    "p3": {"sq0": 0.4, "sq_step": 0.014, "sq_n": 100, "kill": 2.3, "label": 3.4},
    "w1": {"card0": 0.4, "card_step": 0.75, "card_n": 3, "note": 3.3},
    "w2": {"type0": 0.4, "type_rate": 40.0, "note0": 2.0, "note_step": 0.4, "note_n": 3,
           "stage0": 3.9, "stage_step": 0.75, "stage_n": 3, "done": 6.3},
    "w3": {"card0": 0.4, "card_step": 0.55, "card_n": 3, "item_step": 0.12, "verdict": 3.8},
    "o4": {"row0": 0.4, "row_step": 0.42, "row_n": 3, "chip0": 2.0, "chip_step": 0.22, "chip_n": 4},
    "h2": {"bar0": 0.35, "bar_step": 0.026, "bar_n": 44, "pin0": 1.8, "pin_step": 0.13, "pin_n": 6, "count": 3.0},
    "h3": {"row0": 0.35, "row_step": 0.32, "row_n": 5, "note": 2.4},
    "h4": {"sq0": 0.3, "sq_step": 0.22, "sq_n": 6, "fatal": 2.4, "stat": 3.6},
    "g1": {"bar0": 0.35, "bar_step": 0.065, "bar_n": 26, "dot": 2.5, "label": 3.1},
    "g2": {"dots": 0.3, "split": 1.2, "pick": 2.3, "travel": 3.0},
    "g3": {"n0": 0.4, "n_step": 0.038, "n_n": 40, "bar": 2.4, "cand0": 3.0, "cand_step": 0.12},
    "o1": {"type0": 0.3, "type_rate": 36.0, "panel0": 1.8, "panel_step": 0.2, "panel_n": 4},
    "o2": {"type0": 0.3, "type_rate": 44.0, "note0": 1.6, "note_step": 0.32, "note_n": 3},
    "o3": {"card0": 0.3, "card_step": 0.38, "card_n": 3, "item_step": 0.11},
    "end": {"helix": 0.2, "word": 0.75, "rule": 1.2, "pip": 1.7, "url": 2.4},
}

# Terminal lines live here so the typing animation and its keystroke sounds count the same characters.
O1_LINES = ("$ pip install gene-edit-check", "$ gene-edit-check doctor")
O2_LINES = ("$ gene-edit-check run-h5ad \\", "   --h5ad my_cells.h5ad \\",
            "   --gene-symbol PDCD1 \\", "   --target-query \"state == 'exhausted'\"")


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
AUDIT = False
TEXT_BOXES: list = []


def text(c, s, x, y, size, color, *, align="left", alpha=1.0, family=LABEL, bold=True,
         dy=0.0, rot=0.0, track=0.0, fit=True):
    """Draw text with its optical centre on y. `track` adds letter spacing, which the small labels want.

    With `AUDIT` on, every visible run records its box in `TEXT_BOXES` so `audit_overlaps()` can find
    collisions by measurement instead of by eye.
    """
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
    if AUDIT and alpha > 0.25 and rot == 0.0 and total > 0:
        base = y + dy + size * 0.35 * sc
        bx = x + ox * sc
        TEXT_BOXES.append((bx, base - 0.72 * size * sc, bx + total * sc, base + 0.21 * size * sc, s))
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


def arrow(c, x0, y0, x1, y1, w, color, alpha=1.0, head=22):
    """A rule with a solid head, used wherever one thing feeds another."""
    if alpha <= 0:
        return
    ang = math.atan2(y1 - y0, x1 - x0)
    rule(c, x0, y0, x1 - head * 0.85 * math.cos(ang), y1 - head * 0.85 * math.sin(ang), w, color, alpha)
    c.new_path()
    c.set_source_rgba(*color, alpha)
    c.move_to(x1, y1)
    c.line_to(x1 - head * math.cos(ang - 0.42), y1 - head * math.sin(ang - 0.42))
    c.line_to(x1 - head * math.cos(ang + 0.42), y1 - head * math.sin(ang + 0.42))
    c.close_path()
    c.fill()


def leader(c, x0, y0, x1, y1, label, color, alpha=1.0, size=28, align="left", gap=16):
    """A thin elbow from a mark to its name, so nothing on screen goes unexplained."""
    if alpha <= 0:
        return
    rule(c, x0, y0, x1, y0, 3, color, alpha * 0.75)
    rule(c, x1, y0, x1, y1, 3, color, alpha * 0.75)
    text(c, label, x1 + (gap if align == "left" else -gap), y1, size, color, align=align,
         alpha=alpha, family=LABEL, bold=False, fit=False)


# ---- biology ---------------------------------------------------------------------------------
_BASES = "ATCGGCTAACGTTGCAGCTAACGTGGCATCGGATCCGTAAGGCTTCAGGACCTGGATGCAAGGTCAGTACG"
BASE_COLOR = {"A": BLUE, "T": YELLOW, "C": RED, "G": GREY_D}
PAIR = {"A": "T", "T": "A", "C": "G", "G": "C"}


def helix(c, cx, cy, height, amp, rot, *, n=22, turns=2.0, grow=1.0, alpha=1.0, rail=BLACK,
          dot=13.0, seq=_BASES):
    """A turning double helix: two rails, coloured base-pair rungs, depth taken from each rung's phase."""
    items = []
    for i in range(n):
        u = i / max(n - 1, 1)
        if u > grow:
            break
        y = cy - height / 2 + height * u
        ph = u * turns * 2 * math.pi + rot
        items.append((y, cx + amp * math.sin(ph), cx - amp * math.sin(ph), math.cos(ph), seq[i % len(seq)]))
    for y, xa, xb, d, b in items:
        rule(c, xa, y, (xa + xb) / 2, y, 6, BASE_COLOR[b], alpha * (0.45 + 0.55 * (d + 1) / 2))
        rule(c, (xa + xb) / 2, y, xb, y, 6, BASE_COLOR[PAIR[b]], alpha * (0.45 + 0.55 * (-d + 1) / 2))
    for k in range(len(items) - 1):
        y0, xa0, xb0, d0, _ = items[k]
        y1, xa1, xb1, _, _ = items[k + 1]
        rule(c, xa0, y0, xa1, y1, 7, rail, alpha * (0.35 + 0.4 * max(d0, 0)), cairo.LINE_CAP_ROUND)
        rule(c, xb0, y0, xb1, y1, 7, rail, alpha * (0.35 + 0.4 * max(-d0, 0)), cairo.LINE_CAP_ROUND)
    pts = []
    for y, xa, xb, d, b in items:
        pts += [(d, xa, y, BASE_COLOR[b]), (-d, xb, y, BASE_COLOR[PAIR[b]])]
    for d, x, y, col in sorted(pts, key=lambda q: q[0]):
        disc(c, x, y, dot * (0.72 + 0.32 * d), col, alpha * (0.5 + 0.5 * (d + 1) / 2))


def dna_strand(c, x0, x1, y, *, n=24, alpha=1.0, rail=BLACK, reveal=1.0, bases=True, size=20,
               hot=None, hot_letter="T", sep=46, seq=_BASES):
    """An unwound double strand: two rails with base-pair rungs between them."""
    step = (x1 - x0) / max(n - 1, 1)
    end = x0 + (x1 - x0) * clamp(reveal)
    rule(c, x0 - 14, y - sep, end + 14, y - sep, 7, rail, alpha * 0.9)
    rule(c, x0 - 14, y + sep, end + 14, y + sep, 7, rail, alpha * 0.9)
    for i in range(n):
        if i / max(n - 1, 1) > reveal:
            break
        x = x0 + i * step
        is_hot = hot is not None and i == hot
        b = hot_letter if is_hot else seq[i % len(seq)]
        col = RED if is_hot else BASE_COLOR[b]
        rule(c, x, y - sep + 4, x, y + sep - 4, 10 if is_hot else 6, col, alpha * (1.0 if is_hot else 0.85))
        if bases:
            text(c, b, x, y - sep - 26, size, col, align="center", alpha=alpha * (1.0 if is_hot else 0.8),
                 family=MONO, bold=False, fit=False)
            text(c, PAIR[b], x, y + sep + 26, size, BASE_COLOR[PAIR[b]], align="center", alpha=alpha * 0.55,
                 family=MONO, bold=False, fit=False)
    return step


def cell(c, x, y, r, col, *, alpha=1.0, t=0.0, seed=0, level=1.0, dark=True, nucleus=True):
    """A cell: membrane, nucleus, and chromatin marks whose count stands for how active the gene is."""
    if r <= 2 or alpha <= 0:
        return
    disc(c, x, y, r, col, alpha * 0.16)
    ring(c, x, y, r, col, max(4.0, r * 0.055), alpha)
    if not nucleus or r < 30:
        return
    nr, nx, ny = r * 0.36, x - r * 0.16, y - r * 0.08
    disc(c, nx, ny, nr, col, alpha * 0.85)
    rng = np.random.default_rng(seed)
    for _ in range(6):
        keep = rng.uniform(0, 1) <= level
        a = rng.uniform(0, 2 * math.pi) + t * 0.25
        rr = nr * 0.55 * math.sqrt(rng.uniform(0.1, 1))
        if not keep:
            continue
        s = max(nr * 0.17, 5)
        rect(c, nx + rr * math.cos(a) - s / 2, ny + rr * math.sin(a) - s / 2, s, s, PAPER, alpha * 0.9)


def editor(c, x, y, s, col, alpha=1.0, t=0.0):
    """The base editor, drawn as a clamp: two half-discs closing on the base it is working on."""
    if alpha <= 0 or s <= 0:
        return
    open_ = 1.0 + 0.1 * math.sin(t * 3)
    half_disc(c, x, y - s * 0.34 * open_, s, col, math.pi, alpha)
    half_disc(c, x, y + s * 0.34 * open_, s, col, 0.0, alpha)
    rect(c, x - s * 0.1, y - s * 1.6, s * 0.2, s * 0.7, col, alpha)


# ---- chapter furniture ---------------------------------------------------------------------------
def chapter(c, t, num, label, fg, *, t0=0.2):
    a = fade(t, t0, 0.6)
    if a <= 0.004:
        return
    rule(c, W - 430, 112, W - 80, 112, 4, fg, a * 0.9)
    text(c, num, W - 430, 172, 34, fg, alpha=a, family=DISPLAY)
    text(c, label, W - 80, 172, 26, fg, align="right", alpha=a * 0.72, family=LABEL, bold=False, track=5)


def headline(c, t, s, x, y, size, fg, *, t0=0.3, sub=None, sub_gap=76):
    a = fade(t, t0, 0.55)
    dx = (1 - out_cubic((t - t0) / 0.7)) * 36
    text(c, s, x - dx, y, size, fg, alpha=a, family=DISPLAY)
    if sub:
        a2 = fade(t, t0 + 0.4, 0.6)
        text(c, sub, x, y + sub_gap, 32, fg, alpha=a2 * 0.68, family=LABEL, bold=False)


def spine(c, t_abs, paper):
    """A faint strand scrolling under the middle chapters: the film's connective tissue."""
    col = BLACK if paper else PAPER
    y, step, n = 1014, 46.0, 44
    off = (t_abs * 26.0) % step
    rule(c, 0, y - 20, W, y - 20, 4, col, 0.13)
    rule(c, 0, y + 20, W, y + 20, 4, col, 0.13)
    for i in range(n + 2):
        x = -step + i * step + off
        b = _BASES[int(i + t_abs * 26.0 / step) % len(_BASES)]
        rule(c, x, y - 18, x, y + 18, 5, BASE_COLOR[b], 0.2)


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
    return {"cands": cands, "abs": [abs(x["lfc"]) for x in cands], "null": null["abs_lfc_values"],
            "cutoff": null["abs_lfc_max"], "n_edits": len(its), "n_null": null["n_random_edits"],
            "n_fail_meth": sum(1 for it in its if not it["report"]["tier2_hazards"]["methylation_hazard_passed"])}


# ================================================================================================
# TITLE
# ================================================================================================
def scene_title(c, t, dur, paper):
    T = TIMING["title"]
    p = out_cubic(fade(t, T["helix"], 1.0))
    helix(c, 330, 560, 760 * p, 170, 0.4 + t * 0.7, n=20, grow=p, alpha=p, rail=BLACK, dot=15)
    p = out_cubic(fade(t, T["square"], 0.6))
    rect(c, 1470, lerp(-220, 150, p), 190, 190, YELLOW)
    p = out_cubic(fade(t, T["halfdisc"], 0.7))
    half_disc(c, 1330, H + 30, 230 * p, RED, math.pi)
    p = out_cubic(fade(t, T["rule"], 0.7))
    rule(c, 700, 626, lerp(700, 1800, p), 626, 8, BLACK)
    a = fade(t, T["word"], 0.6)
    dx = (1 - out_cubic((t - T["word"]) / 0.8)) * 44
    w = text(c, "gene-edit-", 700 - dx, 540, 104, BLACK, alpha=a, family=DISPLAY)
    text(c, "check", 700 - dx + w, 540, 104, RED, alpha=a, family=DISPLAY)
    a = fade(t, T["sub"], 0.7)
    text(c, "check a gene edit on a computer,", 704, 700, 36, BLACK, alpha=a * 0.8, family=LABEL, bold=False)
    text(c, "before anyone spends weeks on it in a lab.", 704, 750, 36, BLACK, alpha=a * 0.8, family=LABEL, bold=False)
    for k in range(3):
        a2 = out_back(fade(t, T["cells"] + k * 0.22, 0.5), 1.8)
        if a2 <= 0:
            continue
        cell(c, 712 + k * 142, 884 + 10 * math.sin(t * 1.4 + k), 54 * min(a2, 1.0),
             (BLUE, RED, YELLOW)[k], alpha=min(a2, 1.0), t=t, seed=k, dark=False)


# ================================================================================================
# 01  CONTEXT  — why a single-letter edit is worth checking at all
# ================================================================================================
def _c1(c, t):
    """One genome, read differently by every kind of cell."""
    T = TIMING["c1"]
    p = out_cubic(fade(t, T["helix"], 1.1))
    helix(c, 440, 640, 520 * p, 140, 0.2 + t * 0.6, n=18, grow=p, alpha=p, rail=PAPER, dot=13)
    text(c, "one genome", 440, 950, 30, PAPER, align="center", alpha=fade(t, T["helix"] + 0.9) * 0.85,
         family=LABEL, bold=False, track=3)
    names, cols, lv = ["skin", "muscle", "nerve", "T cell"], [YELLOW, RED, BLUE, RED], [0.25, 0.55, 0.9, 1.0]
    for k in range(T["bud_n"]):
        a = min(out_back(fade(t, T["bud0"] + k * T["bud_step"], 0.5), 1.7), 1.0)
        if a <= 0:
            continue
        x, y = 1090 + (k % 2) * 400, 440 + (k // 2) * 290
        arrow(c, 640, 640, x - 118, y, 4, PAPER, a * 0.3, 16)
        cell(c, x, y, 100 * a, cols[k], alpha=a, t=t, seed=k, level=lv[k])
        text(c, names[k], x, y + 142, 28, PAPER, align="center", alpha=a * 0.75,
             family=LABEL, bold=False, track=3)
    a = fade(t, T["read"], 0.7)
    text(c, "same letters everywhere. each cell reads a different part of them.", 1290, 960, 32, PAPER,
         align="center", alpha=a * 0.85, family=LABEL, bold=False)
    headline(c, t, "every cell carries the same DNA", 150, 250, 82, PAPER, t0=0.1,
             sub="what separates a skin cell from a T cell is which parts get read")


def _c2(c, t):
    """A base editor changes one letter, with no cut."""
    T = TIMING["c2"]
    x0, x1, y = 300, 1620, 600
    rev = out_cubic(fade(t, T["strand"], 1.2))
    flipped = t > T["flip"]
    step = dna_strand(c, x0, x1, y, n=22, reveal=rev, rail=PAPER, hot=11 if t > T["editor"] else None,
                      hot_letter="T" if flipped else "C", size=22)
    ed = out_cubic(fade(t, T["editor"], 0.8))
    if ed > 0:
        ex = x0 + 11 * step
        editor(c, ex, y, 84 * ed, YELLOW, ed, t)
        leader(c, ex + 86, y - 150, 1250, 350, "a base editor", YELLOW, fade(t, T["editor"] + 0.4) * 0.95, 30)
    if flipped:
        pr = clamp((t - T["flip"]) / 0.8)
        ring(c, x0 + 11 * step, y, 50 + 130 * pr, RED, 5, 1 - pr)
    a = fade(t, T["label"], 0.6)
    text(c, "C", 640, 880, 56, RED, align="center", alpha=a, family=DISPLAY)
    arrow(c, 690, 880, 810, 880, 6, PAPER, a, 20)
    text(c, "T", 860, 880, 56, YELLOW, align="center", alpha=a, family=DISPLAY)
    text(c, "one base swapped. the strand is never cut.", 960, 880, 32, PAPER,
         alpha=a * 0.85, family=LABEL, bold=False)
    headline(c, t, "a base editor changes one letter", 150, 250, 82, PAPER, t0=0.1,
             sub="enough to turn a gene down without breaking the strand open")


def _c3(c, t):
    """Whatever carries the editor reaches sick and healthy cells alike."""
    T = TIMING["c3"]
    dim = out_cubic((t - T["dim"]) / 1.3)
    a1 = min(out_back(fade(t, T["cell1"], 0.5), 1.6), 1.0)
    a2 = min(out_back(fade(t, T["cell2"], 0.5), 1.6), 1.0)
    cell(c, 620, 650, 200 * a1, RED, alpha=a1, t=t, seed=3, level=lerp(1.0, 0.2, dim))
    cell(c, 1340, 650, 200 * a2, BLUE, alpha=a2, t=t, seed=4, level=lerp(1.0, 0.3, dim))
    text(c, "the cells causing trouble", 620, 900, 30, PAPER, align="center",
         alpha=fade(t, T["cell1"] + 0.4) * 0.85, family=LABEL, bold=False, track=2)
    text(c, "the ones you need to keep", 1340, 900, 30, PAPER, align="center",
         alpha=fade(t, T["cell2"] + 0.4) * 0.85, family=LABEL, bold=False, track=2)
    p = out_cubic(fade(t, T["dim"] - 0.6, 0.9))
    if p > 0:
        arrow(c, 250, 420, lerp(250, 1700, p), 420, 7, YELLOW, 1.0, 24)
        text(c, "the editor goes wherever the carrier goes", 250, 368, 30, YELLOW,
             alpha=p * 0.9, family=LABEL, bold=False)
    ab = fade(t, T["bad"], 0.5)
    if ab > 0:
        s = out_back(ab, 1.9) * 150
        ring(c, 1340, 650, s, RED, 10, ab)
        slash(c, 1340, 650, s, RED, 10, ab)
        text(c, "nobody can see this part in advance", 1340, 952, 30, RED, align="center",
             alpha=ab * 0.95, family=LABEL, bold=False)
    headline(c, t, "an edit cannot pick its cells", 150, 250, 82, PAPER, t0=0.1,
             sub="so a useful edit has to act on one state and spare the other")


def scene_context(c, t, dur, paper):
    chapter(c, t, "01", "WHY", PAPER)
    stage(c, t, [(d, fn) for _, fn, d in BEATS["context"]])


# ================================================================================================
# 02  PROBLEM
# ================================================================================================
def _p1(c, t):
    """Sixty places the editor could act, and no way to tell which is right."""
    T = TIMING["p1"]
    x0, x1, y = 240, 1690, 560
    rev = out_cubic(fade(t, T["strand"], 1.0))
    step = dna_strand(c, x0, x1, y, n=26, reveal=rev, rail=BLACK, size=19, sep=40)
    for k in range(T["pin_n"]):
        p = out_back(fade(t, T["pin0"] + k * T["pin_step"], 0.3), 2.2)
        if p <= 0:
            continue
        x = x0 + [2, 5, 8, 11, 14, 17, 20, 23][k] * step
        rule(c, x, y - 92, x, y - 92 - 88 * p, 4, RED, min(p, 1.0) * 0.8)
        disc(c, x, y - 192 * min(p, 1.0), 15 * min(p, 1.0), RED)
    text(c, "each pin marks a spot the editor can reach", 240, 714, 28, RED,
         alpha=fade(t, T["pin0"] + 0.9) * 0.9, family=LABEL, bold=False)
    a = fade(t, T["count"], 0.6)
    text(c, f"{int(60 * out_cubic((t - T['count']) / 1.1))}", 1690, 850, 140, RED, align="right",
         alpha=a, family=DISPLAY)
    text(c, "places on this one gene", 1690, 934, 30, BLACK, align="right",
         alpha=a * 0.75, family=LABEL, bold=False)
    text(c, "which one turns the gene down without harming anything else?", 240, 900, 34, BLACK,
         alpha=fade(t, T["count"] + 0.5) * 0.9, family=LABEL, bold=False)
    headline(c, t, "sixty possible edits", 150, 250, 92, BLACK, t0=0.1,
             sub="one gene, and every spot a base editor could act on")


def _p2(c, t):
    """Bench time, multiplied by the number of candidates."""
    T = TIMING["p2"]
    y = 470
    for i in range(T["seg_n"]):
        p = out_cubic(fade(t, T["seg0"] + i * T["seg_step"], 0.25))
        if p <= 0:
            continue
        rect(c, 240 + i * 116, y - 32, 98 * p, 64, BLACK if i % 2 else GREY_L)
    weeks = int(clamp((t - T["seg0"]) / (T["seg_step"] * T["seg_n"])) * T["seg_n"])
    text(c, f"{weeks} weeks", 240, y + 112, 44, RED, alpha=fade(t, T["seg0"] + 0.3), family=DISPLAY)
    text(c, "to take a single candidate through a bench", 240, y + 166, 28, BLACK,
         alpha=fade(t, T["seg0"] + 0.6) * 0.7, family=LABEL, bold=False)
    text(c, "x 60", 1810, y + 150, 104, RED, align="right", alpha=fade(t, T["grid"] - 0.35), family=DISPLAY)
    cols, s, gap = 20, 36, 10
    for i in range(T["grid_n"]):
        p = out_back(fade(t, T["grid"] + i * T["grid_step"], 0.22), 2.2)
        if p <= 0:
            continue
        rect(c, 240 + (i % cols) * (s + gap), 740 + (i // cols) * (s + gap), s * min(p, 1.0), s * min(p, 1.0), BLACK)
    text(c, "years of bench work to rule them out one at a time", 1180, 800, 32, BLACK,
         alpha=fade(t, T["grid"] + 1.3) * 0.9, family=LABEL, bold=False)
    headline(c, t, "weeks for every candidate", 150, 250, 92, BLACK, t0=0.1,
             sub="and the list does not get shorter by itself")


def _p3(c, t):
    """The general attrition picture, which is not this project's own result."""
    T = TIMING["p3"]
    cols, s, gap = 20, 38, 13
    gx, gy = 240, 440
    kill = out_cubic((t - T["kill"]) / 1.6)
    for i in range(T["sq_n"]):
        p = out_back(fade(t, T["sq0"] + i * T["sq_step"], 0.2), 2.0)
        if p <= 0:
            continue
        survives = i < 14
        a = 1.0 if survives else lerp(1.0, 0.2, kill)
        sz = s * min(p, 1.0) * (1.0 if survives else lerp(1.0, 0.68, kill))
        rect(c, gx + (i % cols) * (s + gap) + (s - sz) / 2, gy + (i // cols) * (s + gap) + (s - sz) / 2,
             sz, sz, BLUE if survives else GREY_L, a)
    text(c, "100 drugs entering a first trial", gx, gy - 46, 28, BLACK,
         alpha=fade(t, T["sq0"] + 0.5) * 0.75, family=LABEL, bold=False, track=2)
    a = fade(t, T["label"], 0.7)
    text(c, "13.8%", 1300, 500, 118, BLUE, alpha=a, family=DISPLAY)
    text(c, "of them are ever approved", 1300, 582, 30, BLACK, alpha=a * 0.85, family=LABEL, bold=False)
    text(c, "3.4% in oncology", 1300, 658, 44, RED, alpha=fade(t, T["label"] + 0.5), family=DISPLAY)
    text(c, "Wong, Siah & Lo, Biostatistics 2019", 1300, 718, 24, BLACK,
         alpha=fade(t, T["label"] + 0.8) * 0.6, family=LABEL, bold=False)
    text(c, "industry-wide, not this project's data", 1300, 752, 24, BLACK,
         alpha=fade(t, T["label"] + 0.8) * 0.6, family=LABEL, bold=False)
    headline(c, t, "most candidates never make it", 150, 250, 92, BLACK, t0=0.1,
             sub="the cheapest failure is the one found before the bench")


def scene_problem(c, t, dur, paper):
    chapter(c, t, "02", "PROBLEM", BLACK)
    stage(c, t, [(d, fn) for _, fn, d in BEATS["problem"]])


# ================================================================================================
# 03  HOW GENE-EDIT-CHECK WORKS
# ================================================================================================
def _w1(c, t):
    """What a person actually has to supply."""
    T = TIMING["w1"]
    cw, gap = 480, 70
    x0 = (W - (3 * cw + 2 * gap)) / 2
    for i in range(T["card_n"]):
        p = out_cubic(fade(t, T["card0"] + i * T["card_step"], 0.6))
        if p <= 0:
            continue
        x, y = x0 + i * (cw + gap), 400 + (1 - p) * 50
        rect(c, x, y, cw, 68, (RED, YELLOW, BLUE)[i], p)
        text(c, ["single-cell data", "a gene", "two cell labels"][i], x + 26, y + 34, 32,
             PAPER if i != 1 else BLACK, alpha=p, family=DISPLAY)
        rect(c, x, y + 68, cw, 300, PAPER, p * 0.07)
        if i == 0:
            for r in range(4):
                a2 = fade(t, T["card0"] + 0.35 + r * 0.09, 0.3)
                rect(c, x + 30, y + 120 + r * 48, 250, 30, PAPER, a2 * 0.22)
                rect(c, x + 296, y + 120 + r * 48, 150, 30, (RED, BLUE)[r % 2], a2 * 0.75)
            text(c, "one row per cell, one label column", x + 26, y + 338, 22, PAPER,
                 alpha=fade(t, T["card0"] + 0.9) * 0.75, family=LABEL, bold=False, fit=False)
        elif i == 1:
            helix(c, x + cw / 2, y + 196, 160, 44, 0.3 + t * 0.9, n=9,
                  alpha=fade(t, T["card0"] + T["card_step"] + 0.3), rail=PAPER, dot=8)
            text(c, "PDCD1", x + cw / 2, y + 310, 34, YELLOW, align="center",
                 alpha=fade(t, T["card0"] + T["card_step"] + 0.5), family=MONO, bold=False)
        else:
            tb = T["card0"] + 2 * T["card_step"]
            cell(c, x + 150, y + 180, 72, RED, alpha=fade(t, tb + 0.3), t=t, seed=7)
            cell(c, x + 330, y + 180, 72, BLUE, alpha=fade(t, tb + 0.45), t=t, seed=8)
            text(c, "to change", x + 150, y + 290, 26, PAPER, align="center",
                 alpha=fade(t, tb + 0.5) * 0.8, family=LABEL, bold=False, fit=False)
            text(c, "to protect", x + 330, y + 290, 26, PAPER, align="center",
                 alpha=fade(t, tb + 0.6) * 0.8, family=LABEL, bold=False, fit=False)
    text(c, "a public atlas works. so does an experiment of your own.", W / 2, 856, 32, PAPER,
         align="center", alpha=fade(t, T["note"], 0.7) * 0.9, family=LABEL, bold=False)
    text(c, "no reference genome to prepare, no training to run", W / 2, 906, 28, PAPER,
         align="center", alpha=fade(t, T["note"] + 0.4, 0.7) * 0.65, family=LABEL, bold=False)
    headline(c, t, "what you need to start", 150, 250, 86, PAPER, t0=0.1,
             sub="three things, and two of them are just names")


def _w2(c, t):
    """The one command, and what it spends its time on."""
    T = TIMING["w2"]
    rect(c, 170, 400, 1000, 300, PAPER, fade(t, 0.15) * 0.07)
    rule(c, 170, 400, 170, 700, 8, PAPER, fade(t, 0.15))
    shown = int(max(0.0, t - T["type0"]) * T["type_rate"])
    y = 462
    for s_ in O2_LINES:
        n = min(len(s_), shown)
        shown -= max(n, 0)
        if n <= 0:
            break
        text(c, s_[:n], 214, y, 30, PAPER, family=MONO, bold=False, fit=False)
        y += 68
    for i, (lab, col) in enumerate((("your file", YELLOW), ("your gene", RED), ("your two labels", BLUE))):
        p = out_cubic(fade(t, T["note0"] + i * T["note_step"], 0.5))
        if p <= 0:
            continue
        y2 = 530 + i * 68
        rule(c, 1190, y2, lerp(1190, 1260, p), y2, 5, col)
        text(c, lab, 1280, y2, 28, PAPER, alpha=p * 0.9, family=LABEL, bold=False, fit=False)
    stages = ["find the reachable edits", "predict both cell states", "run the six checks"]
    for i in range(T["stage_n"]):
        p = out_cubic(fade(t, T["stage0"] + i * T["stage_step"], 0.55))
        if p <= 0:
            continue
        bw, bx = 500, 180 + i * 540
        rect(c, bx, 790, bw, 22, PAPER, 0.16)
        rect(c, bx, 790, bw * p, 22, (RED, YELLOW, BLUE)[i])
        text(c, stages[i], bx, 852, 26, PAPER, alpha=min(p * 1.4, 1.0) * 0.85, family=LABEL,
             bold=False, fit=False)
    a = fade(t, T["done"], 0.6)
    text(c, "minutes on one laptop and one free API key", W / 2, 940, 30, YELLOW, align="center",
         alpha=a * 0.95, family=LABEL, bold=False)
    headline(c, t, "one command does the run", 150, 250, 86, PAPER, t0=0.1,
             sub="no code to write, and nothing hard-coded to one gene")


def _w3(c, t):
    """What lands on disk, and the decision it supports."""
    T = TIMING["w3"]
    cards = [("a report", ["every candidate", "every check", "the thresholds used"], RED),
             ("the figures", None, BLUE),
             ("a guide to order", ["chr2:241,859,845", "A to G  ·  ABE8e", "protospacer + PAM"], YELLOW)]
    cw, gap = 470, 56
    x0 = (W - (3 * cw + 2 * gap)) / 2
    for i, (title_, rows, col) in enumerate(cards):
        p = out_cubic(fade(t, T["card0"] + i * T["card_step"], 0.5))
        if p <= 0:
            continue
        x, y = x0 + i * (cw + gap), 390 + (1 - p) * 50
        rect(c, x, y, cw, 64, col, p)
        text(c, title_, x + 24, y + 32, 30, PAPER if i != 2 else BLACK, alpha=p, family=DISPLAY)
        rect(c, x, y + 64, cw, 268, PAPER, p * 0.07)
        if rows:
            for j, r in enumerate(rows):
                a2 = fade(t, T["card0"] + i * T["card_step"] + 0.28 + j * T["item_step"], 0.35)
                rect(c, x + 24, y + 118 + j * 66 - 9, 16, 16, col, a2)
                text(c, r, x + 58, y + 118 + j * 66, 26, PAPER, alpha=a2 * 0.85, family=LABEL,
                     bold=False, fit=False)
        else:
            for j in range(3):
                for k in range(4):
                    a2 = fade(t, T["card0"] + i * T["card_step"] + 0.28 + (j * 4 + k) * 0.04, 0.3)
                    rect(c, x + 42 + k * 98, y + 108 + j * 78, 62, 62,
                         RED if (k == 2 and j != 1) else BLUE, a2 * 0.9)
    a = fade(t, T["verdict"], 0.7)
    arrow(c, W / 2, 770, W / 2, 840, 6, YELLOW, a, 22)
    text(c, "a shortlist worth bench time, or a reason there is none", W / 2, 902, 34, YELLOW,
         align="center", alpha=a, family=LABEL, bold=False)
    headline(c, t, "what lands on disk", 150, 250, 86, PAPER, t0=0.1,
             sub="the evidence, not just a verdict you have to trust")


def scene_work(c, t, dur, paper):
    chapter(c, t, "03", "USING IT", PAPER)
    stage(c, t, [(d, fn) for _, fn, d in BEATS["work"]])


_H2_H = [0.14, 0.2, 0.33, 0.25, 0.18, 0.4, 0.62, 0.48, 0.3, 0.22, 0.35, 0.55, 0.78, 0.95, 0.72, 0.5,
         0.36, 0.28, 0.42, 0.6, 0.47, 0.33, 0.25, 0.4, 0.58, 0.44, 0.3, 0.2, 0.28, 0.46, 0.68, 0.52,
         0.35, 0.24, 0.3, 0.44, 0.33, 0.22, 0.28, 0.38, 0.26, 0.18, 0.22, 0.15]
_H2_PINS = [6, 12, 13, 19, 30, 36]


def _h2(c, t):
    """Step one: where the model is sensitive, and which of those spots an editor can reach."""
    T = TIMING["h2"]
    base, x0, bw = 790, 190, 34
    for i in range(T["bar_n"]):
        p = out_cubic(fade(t, T["bar0"] + i * T["bar_step"], 0.3))
        if p <= 0:
            continue
        h = _H2_H[i] * 390 * p
        rect(c, x0 + i * (bw + 4), base - h, bw, h, BLUE if _H2_H[i] >= 0.6 else PAPER,
             1.0 if _H2_H[i] >= 0.6 else 0.5)
    rule(c, x0, base + 12, x0 + T["bar_n"] * (bw + 4) - 4, base + 12, 3, PAPER, fade(t, 0.5) * 0.4)
    for i in range(T["bar_n"]):
        if i / T["bar_n"] > out_cubic((t - 0.55) / 1.3):
            break
        b = _BASES[i % len(_BASES)]
        text(c, b, x0 + i * (bw + 4) + bw / 2, base + 50, 20, BASE_COLOR[b], align="center",
             alpha=0.75, family=MONO, bold=False, fit=False)
    for k in range(T["pin_n"]):
        p = min(out_back(fade(t, T["pin0"] + k * T["pin_step"], 0.3), 2.2), 1.0)
        if p <= 0:
            continue
        i = _H2_PINS[k]
        x = x0 + i * (bw + 4) + bw / 2
        h = _H2_H[i] * 390
        rule(c, x, base - h - 10, x, base - h - 52, 4, RED, p)
        disc(c, x, base - h - 70, 14 * p, RED)
    al = fade(t, 1.0, 0.6)
    rect(c, 190, 892, 30, 26, BLUE, al)
    text(c, "how much the model says this spot matters", 238, 905, 28, PAPER, alpha=al * 0.85,
         family=LABEL, bold=False)
    disc(c, 205, 955, 13, RED, al)
    text(c, "a spot a base editor can actually reach", 238, 955, 28, PAPER, alpha=al * 0.85,
         family=LABEL, bold=False)
    a = fade(t, T["count"], 0.6)
    text(c, f"{int(60 * out_cubic((t - T['count']) / 1.0))}", 1810, 900, 100, RED, align="right",
         alpha=a, family=DISPLAY)
    text(c, "candidates to score", 1810, 978, 28, PAPER, align="right", alpha=a * 0.8,
         family=LABEL, bold=False)
    headline(c, t, "step one  ·  find the reachable edits", 150, 250, 82, PAPER, t0=0.1,
             sub="the model's own gradients say where this gene is sensitive")


def _h3(c, t):
    """Step two: predict every candidate in both cell states, and read the gap."""
    T = TIMING["h3"]
    mid, scale = 1480, 3100.0
    rule(c, mid, 420, mid, 880, 4, PAPER, fade(t, 0.2) * 0.4)
    text(c, "no change", mid, 926, 24, PAPER, align="center", alpha=fade(t, 0.3) * 0.6,
         family=LABEL, bold=False, track=4)
    text(c, "predicted drop in gene activity", mid - 40, 388, 26, PAPER, align="right",
         alpha=fade(t, 0.3) * 0.65, family=LABEL, bold=False)
    for i, cd in enumerate(DATA["cands"][:T["row_n"]]):
        y = 470 + i * 88
        p = out_cubic(fade(t, T["row0"] + i * T["row_step"], 0.55))
        if p <= 0:
            continue
        wt = abs(2 ** cd["lfc"] - 1) * scale * p
        wh = abs(2 ** cd["healthy"] - 1) * scale * p
        rect(c, mid - wt, y - 29, wt, 25, RED)
        rect(c, mid - wh, y + 2, wh, 25, BLUE)
        if wt - wh > 60:
            rect(c, mid + 32, y - 11, 22, 22, YELLOW, fade(t, T["row0"] + i * T["row_step"] + 0.5))
    a = fade(t, T["note"], 0.6)
    rect(c, 200, 492, 32, 25, RED, a)
    text(c, "cells to change", 248, 504, 28, PAPER, alpha=a * 0.85, family=LABEL, bold=False)
    rect(c, 200, 544, 32, 25, BLUE, a)
    text(c, "cells to protect", 248, 556, 28, PAPER, alpha=a * 0.85, family=LABEL, bold=False)
    rect(c, 200, 598, 22, 22, YELLOW, a)
    text(c, "a wide gap: worth keeping", 248, 609, 28, YELLOW, alpha=a * 0.95, family=LABEL, bold=False)
    headline(c, t, "step two  ·  predict both cell states", 150, 250, 82, PAPER, t0=0.1,
             sub="a candidate is only useful when the two bars disagree")


def _h4(c, t):
    """Step three: six checks, any one of which ends a candidate."""
    T = TIMING["h4"]
    names = ["splicing", "methylation", "neighbours", "bystanders", "tissues", "healthy cells"]
    s, gap = 180, 38
    x0 = (W - (T["sq_n"] * s + (T["sq_n"] - 1) * gap)) / 2
    for i in range(T["sq_n"]):
        p = min(out_back(fade(t, T["sq0"] + i * T["sq_step"], 0.3), 2.0), 1.0)
        if p <= 0:
            continue
        bad = i == 1
        x, y = x0 + i * (s + gap), 390
        rect(c, x, y, s, s, PAPER, p * 0.1)
        fill = out_cubic(fade(t, T["sq0"] + i * T["sq_step"] + 0.22, 0.4))
        if fill > 0:
            rect(c, x, y + s * (1 - fill), s, s * fill, RED if bad else BLUE)
        if bad and fill > 0.9:
            slash(c, x + s / 2, y + s / 2, s * 0.3, PAPER, 10, (fill - 0.9) * 10)
        text(c, names[i], x + s / 2, y + s + 42, 26, PAPER, align="center", alpha=p * 0.85,
             family=LABEL, bold=False, fit=False)
    a = fade(t, T["fatal"], 0.6)
    text(c, "any one of the six ends the candidate", W / 2, 700, 42, YELLOW, align="center",
         alpha=a, family=DISPLAY)
    a2 = fade(t, T["stat"], 0.7)
    text(c, "methylation was the fatal one here. at another gene it is usually a different check.",
         W / 2, 782, 30, PAPER, align="center", alpha=a2 * 0.85, family=LABEL, bold=False)
    text(c, f"{DATA['n_fail_meth']} of {DATA['n_edits']} verified candidates died at this step",
         W / 2, 842, 30, RED, align="center", alpha=fade(t, T["stat"] + 0.4) * 0.95, family=LABEL, bold=False)
    headline(c, t, "step three  ·  six checks, all fatal", 150, 250, 82, PAPER, t0=0.1,
             sub="every limit is fixed before a single prediction is made")


def scene_how(c, t, dur, paper):
    chapter(c, t, "04", "INSIDE", PAPER)
    stage(c, t, [(d, fn) for _, fn, d in BEATS["how"]])


# ================================================================================================
# 04  WHAT YOU GET
# ================================================================================================
def _g1(c, t):
    """Bench time against one run."""
    T = TIMING["g1"]
    y = 520
    for i in range(T["bar_n"]):
        p = out_cubic(fade(t, T["bar0"] + i * T["bar_step"], 0.22))
        if p <= 0:
            continue
        rect(c, 240 + i * 61, y - 38, 48 * p, 76 * p, BLACK if i % 2 else GREY_L)
    text(c, "sixty candidates at a bench", 240, y + 104, 30, BLACK, alpha=fade(t, T["bar0"] + 0.4) * 0.8,
         family=LABEL, bold=False, track=2)
    a = min(out_back(fade(t, T["dot"], 0.5), 2.2), 1.0)
    if a > 0:
        rect(c, 240, y + 240, 48 * a, 76 * a, RED)
        text(c, "one run of gene-edit-check", 240, y + 390, 30, RED, alpha=fade(t, T["dot"] + 0.2) * 0.95,
             family=LABEL, bold=False, track=2)
    text(c, "minutes", 1790, y + 280, 110, RED, align="right", alpha=fade(t, T["label"], 0.7), family=DISPLAY)
    headline(c, t, "weeks become minutes", 150, 250, 92, BLACK, t0=0.1,
             sub="the obvious losers never reach a bench at all")


def _g2(c, t):
    """Propose on one half of the patients, judge on the other."""
    T = TIMING["g2"]
    rng = np.random.default_rng(21)
    pts = [(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(44)]
    sp = out_cubic((t - T["split"]) / 1.5)
    for i, (ux, uy) in enumerate(pts):
        a = fade(t, T["dots"] + i * 0.014, 0.4)
        grp = i % 2
        cx = lerp(960, 600 if grp == 0 else 1360, sp)
        col = BLACK if sp < 0.3 else (RED if grp == 0 else BLUE)
        rect(c, cx + ux * lerp(320, 185, sp) - 11, 630 + uy * lerp(230, 185, sp) - 11, 22, 22, col, a)
    if sp > 0.05:
        rule(c, 960, 390, 960, 930, 6, YELLOW, sp)
        text(c, "split by patient, before anything runs", 960, 976, 26, YELLOW, align="center",
             alpha=sp, family=LABEL, bold=False, track=3)
    text(c, "edits proposed here", 600, 912, 28, RED, align="center", alpha=sp * 0.95, family=LABEL, bold=False)
    text(c, "and graded here", 1360, 912, 28, BLUE, align="center", alpha=sp * 0.95, family=LABEL, bold=False)
    pick = min(out_back(fade(t, T["pick"], 0.4), 2.0), 1.0)
    if pick > 0:
        rect(c, 600 - 26 * pick, 630 - 26 * pick, 52 * pick, 52 * pick, YELLOW)
    tr = out_cubic((t - T["travel"]) / 1.3)
    if tr > 0:
        x = lerp(600, 1360, tr)
        arrow(c, 600, 630, x, 630, 4, YELLOW, 0.7, 18)
        rect(c, x - 26, 604, 52, 52, YELLOW)
        if tr > 0.95:
            ring(c, 1360, 630, 185, YELLOW, 5, clamp((tr - 0.95) * 20))
    headline(c, t, "judged on patients held back", 150, 250, 92, BLACK, t0=0.1,
             sub="the half that proposed an edit never gets to grade it")


def _g3(c, t):
    """Nothing counts until it clears a crowd of random edits."""
    T = TIMING["g3"]
    ax0, ax1, ay = 240, 1700, 790
    rule(c, ax0, ay, ax1, ay, 4, BLACK, fade(t, 0.3) * 0.5)
    xmax = max(DATA["abs"]) * 1.12
    jit = np.random.default_rng(77).normal(0, 13, len(DATA["null"]))

    def px(v):
        return ax0 + 40 + math.sqrt(min(v, xmax) / xmax) * (ax1 - ax0 - 80)

    for i, v in enumerate(DATA["null"]):
        p = min(out_back(fade(t, T["n0"] + i * T["n_step"], 0.25), 2.2), 1.0)
        if p <= 0:
            continue
        rect(c, px(v) + jit[i] - 11 * p, ay - 40 - 34 * (i % 5) - 11 * p, 22 * p, 22 * p, GREY_L)
    text(c, f"{DATA['n_null']} random edits", ax0, 640, 30, GREY_D,
         alpha=fade(t, T["n0"] + 0.5) * 0.95, family=LABEL, bold=False, track=2)
    text(c, "scored at the same gene, first", ax0, 686, 26, GREY_D,
         alpha=fade(t, T["n0"] + 0.8) * 0.8, family=LABEL, bold=False)
    a = fade(t, T["bar"], 0.6)
    cx = px(DATA["cutoff"])
    rule(c, cx, ay + 40, cx, ay - 280, 6, YELLOW, a)
    text(c, "the bar: the best of them", cx + 24, ay - 308, 28, YELLOW, alpha=a, family=LABEL, bold=False)
    for i, v in enumerate(DATA["abs"]):
        p = min(out_back(fade(t, T["cand0"] + i * T["cand_step"], 0.3), 2.0), 1.0)
        if p <= 0:
            continue
        s = 34 * p
        rect(c, px(v) - s / 2, ay - 60 - 50 * (i % 3) - s / 2, s, s, RED if v > DATA["cutoff"] else GREY_D)
    text(c, "only a candidate past the bar is reported at all", ax1, 886, 30, BLACK, align="right",
         alpha=fade(t, T["cand0"] + 1.0) * 0.9, family=LABEL, bold=False)
    headline(c, t, "every claim has to beat chance", 150, 250, 92, BLACK, t0=0.1,
             sub="forty random edits set the bar before any candidate is scored")


def scene_gain(c, t, dur, paper):
    chapter(c, t, "05", "WHY TRUST IT", BLACK)
    stage(c, t, [(d, fn) for _, fn, d in BEATS["gain"]])


# ================================================================================================
# 05  GET IT
# ================================================================================================
def _o1(c, t):
    T = TIMING["o1"]
    half_disc(c, 1700, H + 30, 280 * out_cubic(fade(t, 0.3, 0.9)), BLUE, math.pi)
    rect(c, 170, 420, 900, 280, PAPER, fade(t, 0.15) * 0.07)
    rule(c, 170, 420, 170, 700, 8, PAPER, fade(t, 0.15))
    shown = int(max(0.0, t - T["type0"]) * T["type_rate"])
    y = 498
    for s in O1_LINES:
        n = min(len(s), shown)
        shown -= max(n, 0)
        if n <= 0:
            break
        text(c, s[:n], 218, y, 36, PAPER, family=MONO, bold=False, fit=False)
        y += 78
    for i, (k, col) in enumerate((("api key", RED), ("reference data", BLUE), ("cache", BLUE), ("ready", YELLOW))):
        p = min(out_back(fade(t, T["panel0"] + i * T["panel_step"], 0.3), 2.0), 1.0)
        if p <= 0:
            continue
        y2 = 430 + i * 82
        rect(c, 1240, y2, 50 * p, 50 * p, col)
        text(c, k, 1316, y2 + 25, 30, PAPER, alpha=p * 0.85, family=LABEL, bold=False)
    headline(c, t, "install gene-edit-check", 150, 250, 86, PAPER, t0=0.1,
             sub="one command, then it reports whatever is missing")


def _o4(c, t):
    """Who the tool is for, and what it costs to try."""
    T = TIMING["o4"]
    rows = [("computational biologists", "with a single-cell dataset already in hand", RED),
            ("gene-editing teams", "deciding which candidate is worth a bench", YELLOW),
            ("anyone screening a locus", "who needs the reasoning written down", BLUE)]
    for i, (who, why, col) in enumerate(rows):
        p = out_cubic(fade(t, T["row0"] + i * T["row_step"], 0.6))
        if p <= 0:
            continue
        y = 420 + i * 120 + (1 - p) * 30
        rect(c, 180, y - 26, 52, 52, col, p)
        text(c, who, 266, y - 10, 36, PAPER, alpha=p, family=DISPLAY)
        text(c, why, 266, y + 36, 27, PAPER, alpha=p * 0.7, family=LABEL, bold=False)
    for i, chip in enumerate(("macOS", "Linux", "Windows", "MIT licensed")):
        p = min(out_back(fade(t, T["chip0"] + i * T["chip_step"], 0.4), 1.8), 1.0)
        if p <= 0:
            continue
        x = 1240 + (i % 2) * 300
        y = 430 + (i // 2) * 120
        rect(c, x, y, 260 * p, 78 * p, PAPER, p * 0.9, stroke=3)
        text(c, chip, x + 130, y + 39, 28, PAPER, align="center", alpha=p, fit=False)
    text(c, "free to run. the model needs a free non-commercial key.", 1240, 760, 28, PAPER,
         alpha=fade(t, T["chip0"] + 1.0) * 0.75, family=LABEL, bold=False)
    headline(c, t, "who this is for", 150, 250, 86, PAPER, t0=0.1,
             sub="anyone who has to choose which edit deserves a bench")


def scene_outro(c, t, dur, paper):
    chapter(c, t, "06", "GET IT", PAPER)
    stage(c, t, [(d, fn) for _, fn, d in BEATS["outro"]])


def scene_end(c, t, dur, paper):
    T = TIMING["end"]
    p = out_cubic(fade(t, T["helix"], 1.0))
    helix(c, 1560, 520, 780 * p, 150, 0.3 + t * 0.7, n=18, grow=p, alpha=p, rail=BLACK, dot=14)
    half_disc(c, 420, H, 200 * out_cubic(fade(t, T["helix"] + 0.3, 0.7)), YELLOW, math.pi)
    a = fade(t, T["word"], 0.6)
    dx = (1 - out_cubic((t - T["word"]) / 0.8)) * 40
    w = text(c, "gene-edit-", 150 - dx, 450, 104, BLACK, alpha=a, family=DISPLAY)
    text(c, "check", 150 - dx + w, 450, 104, RED, alpha=a, family=DISPLAY)
    rule(c, 150, 540, lerp(150, 1120, out_cubic(fade(t, T["rule"], 0.7))), 540, 8, BLACK)
    text(c, "pip install gene-edit-check", 150, 632, 44, RED, alpha=fade(t, T["pip"], 0.6),
         family=MONO, bold=False)
    a = fade(t, T["url"], 0.6)
    text(c, "github.com/aravmendiratta/gene-edit-check", 150, 722, 32, BLACK, alpha=a * 0.8,
         family=LABEL, bold=False)
    for i, col in enumerate((RED, BLUE, YELLOW)):
        rect(c, 150 + i * 56, 810, 40, 40, col, a)


# (key, draw function, seconds). A scene's length above must equal the sum of its beats.
BEATS = {
    "context": [("c1", _c1, 5.9), ("c2", _c2, 5.9), ("c3", _c3, 5.98)],
    "problem": [("p1", _p1, 5.18), ("p2", _p2, 5.18), ("p3", _p3, 5.2)],
    "work": [("w1", _w1, 7.4), ("w2", _w2, 7.4), ("w3", _w3, 7.42)],
    "how": [("h2", _h2, 5.9), ("h3", _h3, 5.9), ("h4", _h4, 5.98)],
    "gain": [("g1", _g1, 5.9), ("g2", _g2, 5.9), ("g3", _g3, 5.98)],
    "outro": [("o1", _o1, 4.44), ("o4", _o4, 4.45)],
}


def beat_at(scene: str, key: str) -> float:
    """Seconds from the start of `scene` to the start of beat `key`."""
    off = 0.0
    for k, _, d in BEATS[scene]:
        if k == key:
            return off
        off += d
    raise KeyError(f"{scene}/{key}")


SCENE_FN = {"title": scene_title, "context": scene_context, "problem": scene_problem,
            "work": scene_work, "how": scene_how, "gain": scene_gain, "outro": scene_outro,
            "end": scene_end}


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


def draw_scene(c, idx, t, t_abs):
    name, bars, paper = SCENES[idx]
    c.set_source_rgb(*(PAPER if paper else BLACK))
    c.paint()
    if name in SPINE_SCENES:
        spine(c, t_abs, paper)
    SCENE_FN[name](c, t, bars * BAR, paper)


TRANS_DUR = 0.75


def render(t, surf):
    c = cairo.Context(surf)
    c.scale(OW / W, OH / H)
    idx = max(i for i in range(len(SCENES)) if STARTS[i] <= t + 1e-9)
    local = t - STARTS[idx]
    if idx > 0 and local < TRANS_DUR:
        draw_scene(c, idx - 1, SCENES[idx - 1][1] * BAR + local, t)
        c.save()
        clip_shape(c, TRANSITIONS[(idx - 1) % len(TRANSITIONS)], local / TRANS_DUR)
        draw_scene(c, idx, local, t)
        c.restore()
    else:
        draw_scene(c, idx, local, t)


# ================================================================================================
# sound
# ================================================================================================
SR = 44100
SCALE = [220.00, 261.63, 293.66, 329.63, 392.00, 440.00, 523.25, 587.33, 659.25]  # A minor


def build_events() -> list[tuple[float, str, float]]:
    """Every sound effect, mirroring the animation schedules in TIMING and the beat table in BEATS."""
    ev: list[tuple[float, str, float]] = []
    S = {name: float(STARTS[i]) for i, (name, _, _) in enumerate(SCENES)}

    def add(scene, t, kind, pitch=0.0):
        ev.append((S[scene] + t, kind, pitch))

    T = TIMING["title"]
    add("title", T["helix"], "sweep")
    add("title", T["square"], "thunk", SCALE[1])
    add("title", T["halfdisc"], "thunk", SCALE[0])
    add("title", T["rule"], "whoosh")
    add("title", T["word"], "chime", SCALE[5])
    add("title", T["sub"], "bop", SCALE[2])
    for k in range(3):
        add("title", T["cells"] + k * 0.22, "bop", SCALE[k + 3])

    b, T = beat_at("context", "c1"), TIMING["c1"]
    add("context", b + T["helix"], "sweep")
    for k in range(T["bud_n"]):
        add("context", b + T["bud0"] + k * T["bud_step"], "thunk", SCALE[k % 5])
    add("context", b + T["read"], "bop", SCALE[4])
    b, T = beat_at("context", "c2"), TIMING["c2"]
    add("context", b + T["strand"], "whoosh")
    add("context", b + T["editor"], "thunk", SCALE[3])
    add("context", b + T["flip"], "beep", SCALE[6])
    add("context", b + T["label"], "bop", SCALE[2])
    b, T = beat_at("context", "c3"), TIMING["c3"]
    add("context", b + T["cell1"], "thunk", SCALE[2])
    add("context", b + T["cell2"], "thunk", SCALE[0])
    add("context", b + T["dim"] - 0.6, "whoosh")
    add("context", b + T["dim"], "fall", SCALE[1])
    add("context", b + T["bad"], "alarm", SCALE[1])

    b, T = beat_at("problem", "p1"), TIMING["p1"]
    add("problem", b + T["strand"], "whoosh")
    for k in range(T["pin_n"]):
        add("problem", b + T["pin0"] + k * T["pin_step"], "beep", SCALE[k % 6] * 2)
    add("problem", b + T["count"], "chime", SCALE[4])
    b, T = beat_at("problem", "p2"), TIMING["p2"]
    for i in range(T["seg_n"]):
        add("problem", b + T["seg0"] + i * T["seg_step"], "bop", SCALE[i % 5] / 2)
    for i in range(0, T["grid_n"], 4):
        add("problem", b + T["grid"] + i * T["grid_step"], "tick", SCALE[(i // 4) % 6] * 2)
    b, T = beat_at("problem", "p3"), TIMING["p3"]
    for i in range(0, T["sq_n"], 7):
        add("problem", b + T["sq0"] + i * T["sq_step"], "tick", SCALE[(i // 7) % 7] * 2)
    add("problem", b + T["kill"], "fall", SCALE[0])
    add("problem", b + T["label"], "chime", SCALE[3])

    b, T = beat_at("work", "w1"), TIMING["w1"]
    for i in range(T["card_n"]):
        add("work", b + T["card0"] + i * T["card_step"], "thunk", SCALE[i % 5])
        for j in range(3):
            add("work", b + T["card0"] + i * T["card_step"] + 0.35 + j * 0.09, "tick", SCALE[(i + j) % 7] * 2)
    add("work", b + T["note"], "bop", SCALE[4])
    b, T = beat_at("work", "w2"), TIMING["w2"]
    for i in range(int(sum(len(x) for x in O2_LINES) / T["type_rate"] / 0.1)):
        add("work", b + T["type0"] + i * 0.1, "key", SCALE[i % 5] * 2)
    for i in range(T["note_n"]):
        add("work", b + T["note0"] + i * T["note_step"], "swipe", SCALE[i % 5])
    for i in range(T["stage_n"]):
        add("work", b + T["stage0"] + i * T["stage_step"], "fill", SCALE[i % 6] * 2)
        add("work", b + T["stage0"] + i * T["stage_step"] + 0.5, "beep", SCALE[(i + 3) % 6] * 2)
    add("work", b + T["done"], "chime", SCALE[5])
    b, T = beat_at("work", "w3"), TIMING["w3"]
    for i in range(T["card_n"]):
        add("work", b + T["card0"] + i * T["card_step"], "thunk", SCALE[i % 5])
        for j in range(3):
            add("work", b + T["card0"] + i * T["card_step"] + 0.28 + j * T["item_step"], "tick", SCALE[(i + j) % 7] * 2)
    add("work", b + T["verdict"], "chime", SCALE[6])

    b, T = beat_at("how", "h2"), TIMING["h2"]
    for i in range(0, T["bar_n"], 2):
        add("how", b + T["bar0"] + i * T["bar_step"], "tick", SCALE[i % 7] * 2)
    for k in range(T["pin_n"]):
        add("how", b + T["pin0"] + k * T["pin_step"], "beep", SCALE[k % 6] * 2)
    add("how", b + T["count"], "chime", SCALE[5])
    b, T = beat_at("how", "h3"), TIMING["h3"]
    for i in range(T["row_n"]):
        add("how", b + T["row0"] + i * T["row_step"], "swipe", SCALE[i % 5])
    add("how", b + T["note"], "bop", SCALE[3])
    b, T = beat_at("how", "h4"), TIMING["h4"]
    for i in range(T["sq_n"]):
        add("how", b + T["sq0"] + i * T["sq_step"], "bop", SCALE[i % 6])
        add("how", b + T["sq0"] + i * T["sq_step"] + 0.22, "fill", SCALE[i % 6] * 2)
    add("how", b + T["fatal"], "alarm", SCALE[0])
    add("how", b + T["stat"], "bop", SCALE[4])

    b, T = beat_at("gain", "g1"), TIMING["g1"]
    for i in range(0, T["bar_n"], 2):
        add("gain", b + T["bar0"] + i * T["bar_step"], "tick", SCALE[i % 6] * 2)
    add("gain", b + T["dot"], "thunk", SCALE[4])
    add("gain", b + T["label"], "chime", SCALE[6])
    b, T = beat_at("gain", "g2"), TIMING["g2"]
    for i in range(0, 44, 4):
        add("gain", b + T["dots"] + i * 0.014, "tick", SCALE[(i // 4) % 7] * 2)
    add("gain", b + T["split"], "whoosh")
    add("gain", b + T["pick"], "beep", SCALE[5])
    add("gain", b + T["travel"], "swipe", SCALE[3])
    add("gain", b + T["travel"] + 1.25, "chime", SCALE[6])
    b, T = beat_at("gain", "g3"), TIMING["g3"]
    for i in range(0, T["n_n"], 3):
        add("gain", b + T["n0"] + i * T["n_step"], "tick", SCALE[(i // 3) % 7] * 2)
    add("gain", b + T["bar"], "thunk", SCALE[4])
    for i in range(8):
        add("gain", b + T["cand0"] + i * T["cand_step"], "bop", SCALE[i % 6])

    b, T = beat_at("outro", "o1"), TIMING["o1"]
    for i in range(int(sum(len(x) for x in O1_LINES) / T["type_rate"] / 0.1)):
        add("outro", b + T["type0"] + i * 0.1, "key", SCALE[i % 5] * 2)
    for i in range(T["panel_n"]):
        add("outro", b + T["panel0"] + i * T["panel_step"], "bop", SCALE[i % 5])
    b, T = beat_at("outro", "o4"), TIMING["o4"]
    for i in range(T["row_n"]):
        add("outro", b + T["row0"] + i * T["row_step"], "bop", SCALE[i % 5])
    for i in range(T["chip_n"]):
        add("outro", b + T["chip0"] + i * T["chip_step"], "tick", SCALE[i % 6] * 2)

    T = TIMING["end"]
    add("end", T["helix"], "sweep")
    add("end", T["word"], "chime", SCALE[5])
    add("end", T["rule"], "whoosh")
    add("end", T["pip"], "key", SCALE[6])
    add("end", T["url"], "bop", SCALE[4])

    for i in range(1, len(SCENES)):
        t0 = float(STARTS[i])
        ev.append((t0 - TRANS_DUR, "sweep", 0.0))
        ev.append((t0, "impact", 0.0))
    return [e for e in ev if 0 <= e[0] < TOTAL]


def audit_overlaps(step: float = 0.12, min_frac: float = 0.10) -> list[tuple[float, str, str, float]]:
    """Sample the film and report pairs of text runs whose boxes overlap, so collisions are measured.

    Samples skip the crossfades between beats and the transitions between chapters, where two scenes are
    deliberately on screen at once.
    """
    global AUDIT
    surf = cairo.ImageSurface(cairo.FORMAT_ARGB32, OW, OH)
    beat_edges = []
    for i, (name, _, _) in enumerate(SCENES):
        for key, _, _ in BEATS.get(name, []):
            beat_edges.append(float(STARTS[i]) + beat_at(name, key))
    hits: list[tuple[float, str, str, float]] = []
    AUDIT = True
    try:
        t = 0.0
        while t < TOTAL:
            near_scene = any(abs(t - float(x)) < TRANS_DUR + 0.1 for x in STARTS)
            near_beat = any(abs(t - e) < 0.62 for e in beat_edges)
            if near_scene or near_beat:
                t += step
                continue
            TEXT_BOXES.clear()
            render(t, surf)
            boxes = list(TEXT_BOXES)
            for i in range(len(boxes)):
                for j in range(i + 1, len(boxes)):
                    ax0, ay0, ax1, ay1, as_ = boxes[i]
                    bx0, by0, bx1, by1, bs = boxes[j]
                    ow = min(ax1, bx1) - max(ax0, bx0)
                    oh = min(ay1, by1) - max(ay0, by0)
                    if ow <= 0 or oh <= 0:
                        continue
                    small = min((ax1 - ax0) * (ay1 - ay0), (bx1 - bx0) * (by1 - by0))
                    frac = (ow * oh) / max(small, 1.0)
                    if frac >= min_frac:
                        hits.append((round(t, 2), as_, bs, round(frac, 2)))
            t += step
    finally:
        AUDIT = False
        TEXT_BOXES.clear()
    return hits


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
