"""A calm, 41-second motion-graphics explainer for the tool, drawn entirely in code.

    pip install pycairo imageio-ffmpeg numpy scipy
    PYTHONPATH=src python scripts/make_launch_video.py        # writes assets/launch.mp4 (1280x720)

No screenshots and no stock footage. The four steps are each explained with their own animation. The numbers
(the random-edit null, the candidates, the patient and cell counts) are read from the saved results in
`outputs/`. The music is synthesised here (100 BPM, A minor, soft pads and bells), so nothing needs licensing.
Picture and sound are muxed in one ffmpeg pass with the index at the front, so the file streams without stalling.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import cairo
import numpy as np

W, H = 1920, 1080          # the layout grid
OW, OH, FPS = 1280, 720, 30  # what is rendered
BPM = 100
BEAT = 60 / BPM
BAR = BEAT * 4
OUT = Path("outputs")

# A quieter palette than the first cut: deep forest, warm paper, soft teal, one muted coral.
FOREST = (0.055, 0.130, 0.115)
PAPER = (0.957, 0.933, 0.875)
INK = (0.090, 0.150, 0.130)
TEAL = (0.180, 0.620, 0.530)
MINT = (0.610, 0.840, 0.750)
CORAL = (0.929, 0.450, 0.310)
GOLD = (0.910, 0.770, 0.420)
MUTE_D = (0.600, 0.690, 0.650)   # secondary text on dark
MUTE_L = (0.380, 0.450, 0.420)   # secondary text on paper
CARD = (0.985, 0.972, 0.935)

SCENES = [  # name, bars, dark background?
    ("hook", 2, True), ("find", 2, False), ("predict", 2, True), ("stress", 2, False),
    ("judge", 3, True), ("stats", 2, False), ("byo", 2, True), ("end", 2, True),
]
STARTS = np.cumsum([0] + [b * BAR for _, b, _ in SCENES])
TOTAL = float(STARTS[-1])
DATA: dict = {}
SERIF, SANS, MONO = "Georgia", "Segoe UI", "Consolas"


# ---- easing and helpers ----------------------------------------------------------------------
def clamp(x: float) -> float:
    return max(0.0, min(1.0, x))


def out_cubic(x: float) -> float:
    return 1 - (1 - clamp(x)) ** 3


def in_out(x: float) -> float:
    x = clamp(x)
    return x * x * (3 - 2 * x)


def lerp(a: float, b: float, t: float) -> float:
    return a + (b - a) * t


def fade(t: float, t0: float, dur: float = 0.9) -> float:
    return in_out((t - t0) / dur)


def rise(t: float, t0: float, dur: float = 1.0, dist: float = 28.0) -> tuple[float, float]:
    p = out_cubic((t - t0) / dur)
    return p, (1 - p) * dist


def text(c, s, x, y, size, color, *, align="center", alpha=1.0, family=SANS, bold=True, dy=0.0, rot=0.0, fit=True):
    if alpha <= 0.003:
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
    c.set_source_rgba(*color, alpha)
    c.arc(x, y, r, 0, 2 * math.pi)
    c.fill()


def ring(c, x, y, r, color, w, alpha=1.0):
    if r <= 0 or alpha <= 0:
        return
    c.set_source_rgba(*color, alpha)
    c.set_line_width(w)
    c.arc(x, y, r, 0, 2 * math.pi)
    c.stroke()


def line(c, x0, y0, x1, y1, w, color, alpha=1.0):
    c.set_source_rgba(*color, alpha)
    c.set_line_width(w)
    c.set_line_cap(cairo.LINE_CAP_ROUND)
    c.move_to(x0, y0)
    c.line_to(x1, y1)
    c.stroke()


def rrect(c, x, y, w, h, r, color, alpha=1.0, stroke=0.0):
    c.new_sub_path()
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
    c.set_source_rgba(*color, clamp(p * 1.5))
    c.set_line_width(w)
    c.set_line_cap(cairo.LINE_CAP_ROUND)
    c.set_line_join(cairo.LINE_JOIN_ROUND)
    c.move_to(x - 22 * p, y + 2 * p)
    c.line_to(x - 6 * p, y + 18 * p)
    c.line_to(x + 24 * p, y - 18 * p)
    c.stroke()


def cross(c, x, y, p, color, w=9):
    c.set_source_rgba(*color, clamp(p * 1.5))
    c.set_line_width(w)
    c.set_line_cap(cairo.LINE_CAP_ROUND)
    c.move_to(x - 18 * p, y - 18 * p)
    c.line_to(x + 18 * p, y + 18 * p)
    c.move_to(x + 18 * p, y - 18 * p)
    c.line_to(x - 18 * p, y + 18 * p)
    c.stroke()


def base_color(ch: str, dark: bool):
    return {"A": MINT if dark else TEAL, "T": GOLD if dark else (0.78, 0.60, 0.20), "C": CORAL, "G": PAPER if dark else INK}[ch]


PAIR = {"A": "T", "T": "A", "C": "G", "G": "C"}


def helix(c, cx, cy, height, amp, rot, *, alpha=1.0, n=26, turns=2.2, grow=1.0, seq="ATCGGCTAACGTTGCAGCTA", dot=11.0, dark=True):
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


def counter(t, t0, dur, to):
    return to * out_cubic((t - t0) / dur)


def header(c, t, n, title, dark, sub=None):
    fg = PAPER if dark else INK
    mute = MUTE_D if dark else MUTE_L
    a, dy = rise(t, 0.2, 1.0)
    disc(c, 150, 150 + dy, 34, CORAL, a)
    text(c, str(n), 150, 150 + dy, 38, PAPER, alpha=a)
    text(c, title, 215, 150, 78, fg, align="left", alpha=a, family=SERIF, bold=False, dy=dy)
    if sub:
        a2, dy2 = rise(t, 1.0, 1.0)
        text(c, sub, 215, 232, 36, mute, align="left", alpha=a2, bold=False, dy=dy2)


# ---- data ------------------------------------------------------------------------------------
def load_data() -> dict:
    rep = json.loads((OUT / "report_crc-pd1-exhaustion.json").read_text(encoding="utf-8"))
    its = rep["iterations"]
    null = rep["efficacy_null"]
    return {
        "cands": [abs(it["report"]["tier1_selectivity"]["target_delta_sketch_b"]) for it in its],
        "null": null["abs_lfc_values"], "cutoff": null["abs_lfc_max"], "n_edits": len(its),
        "n_fail_meth": sum(1 for it in its if not it["report"]["tier2_hazards"]["methylation_hazard_passed"]),
    }


# ---- scenes (t = seconds since the scene started) --------------------------------------------
def scene_hook(c, t, dark):
    helix(c, 1420, 560, 1000, 230, 0.5 + t * 0.55, grow=out_cubic(t / 2.4), dot=15, n=28)
    a, dy = rise(t, 0.4, 1.4)
    text(c, "One letter of DNA", 140, 450, 112, PAPER, align="left", alpha=a, family=SERIF, bold=False, dy=dy)
    a, dy = rise(t, 1.9, 1.4)
    text(c, "can turn a gene down.", 140, 590, 112, MINT, align="left", alpha=a, family=SERIF, bold=False, dy=dy)
    a, dy = rise(t, 3.5, 1.0)
    text(c, "Here is how we check it before the lab.", 144, 730, 40, MUTE_D, align="left", alpha=a, bold=False, dy=dy)


def smooth_curve(n=240):
    xs = np.linspace(0, 1, n)
    bumps = [(0.50, 0.95, 0.045), (0.43, 0.55, 0.03), (0.60, 0.62, 0.035), (0.29, 0.35, 0.03), (0.74, 0.4, 0.04), (0.14, 0.25, 0.03), (0.88, 0.2, 0.03)]
    y = np.full(n, 0.06)
    for mu, h, s in bumps:
        y += h * np.exp(-((xs - mu) ** 2) / (2 * s * s))
    return xs, y / y.max()


def scene_find(c, t, dark):
    header(c, t, 1, "Find where an edit could work", dark, "The model's gradients show where the gene is most sensitive.")
    x0, x1, base = 160, 1760, 880
    seq = "ACGTTGGCATCGGATCCGTAAGGCTTCAGGACCTGGATGCAAGGTCAGTACGGTAC" * 2
    n = 84
    for i in range(n):
        a = fade(t, 0.4 + i * 0.012, 0.5)
        col = base_color(seq[i], dark)
        x = x0 + i * (x1 - x0) / (n - 1)
        rrect(c, x - 6, base - 22 - 6 * math.sin(i * 0.5), 12, 44, 5, col, a * 0.9)
    xs, ys = smooth_curve()
    prog = out_cubic((t - 1.0) / 2.2)
    top, hmax = base - 70, 380
    pts = [(x0 + u * (x1 - x0), top - v * hmax) for u, v in zip(xs, ys)]
    cut = x0 + prog * (x1 - x0)
    c.save()
    c.rectangle(0, 0, cut, H)
    c.clip()
    c.move_to(pts[0][0], top)
    for p in pts:
        c.line_to(*p)
    c.line_to(pts[-1][0], top)
    c.close_path()
    c.set_source_rgba(*TEAL, 0.35)
    c.fill()
    c.move_to(*pts[0])
    for p in pts[1:]:
        c.line_to(*p)
    c.set_source_rgba(*TEAL, 1)
    c.set_line_width(5)
    c.set_line_join(cairo.LINE_JOIN_ROUND)
    c.stroke()
    c.restore()
    peaks = [0.14, 0.29, 0.43, 0.50, 0.60, 0.74, 0.88]
    for k, u in enumerate(peaks):
        tp = 1.0 + 2.2 * u + 0.15
        a, dy = rise(t, tp, 0.7, 30)
        i = min(int(u * (len(xs) - 1)), len(xs) - 1)
        px, py = pts[i]
        line(c, px, py - 8, px, py - 54 + dy * 0, 3, INK, a * 0.55)
        disc(c, px, py - 62 + dy, 13, CORAL, a)
    a, dy = rise(t, 3.3, 1.0)
    text(c, "each pin: an edit the base editor can reach", 160, 960, 34, MUTE_L, align="left", alpha=a, bold=False, dy=dy)
    v = counter(t, 3.2, 1.2, 60)
    text(c, f"{v:.0f} candidates", 1760, 960, 44, INK, align="right", alpha=fade(t, 3.2), family=SERIF, bold=False)


def cell(c, x, y, r, col, dark, alpha=1.0):
    disc(c, x, y, r, col, alpha * 0.14)
    ring(c, x, y, r, col, 5, alpha)
    disc(c, x - r * 0.18, y - r * 0.12, r * 0.30, col, alpha * 0.35)


def scene_predict(c, t, dark):
    header(c, t, 2, "Predict what each edit does", dark, "AlphaGenome 2 reads the DNA and predicts the effect in each cell type.")
    # the edited letter
    flip = in_out((t - 1.6) / 0.5)
    ch = "C" if flip < 0.5 else "T"
    a = fade(t, 0.5)
    sy = abs(math.cos(flip * math.pi)) if 0 < flip < 1 else 1.0
    c.save()
    c.translate(330, 620)
    c.scale(1, max(sy, 0.03))
    text(c, ch, 0, 0, 220, CORAL, alpha=a, family=SERIF, bold=False)
    c.restore()
    text(c, "one base, edited", 330, 800, 34, MUTE_D, alpha=a, bold=False)
    pr = (t - 2.0) / 1.6
    for k, (cx, col, lab, final) in enumerate(((1130, CORAL, "target cells", 0.35), (1560, MINT, "healthy cells", 1.0))):
        ca = fade(t, 0.8 + k * 0.2)
        cell(c, cx, 600, 150, col, dark, ca)
        rng = np.random.default_rng(3 + k)
        for j in range(26):
            ang, rr = rng.uniform(0, 2 * math.pi), 118 * math.sqrt(rng.uniform(0.15, 1))
            keep = rng.uniform(0, 1) < final
            fade_out = 1.0 if keep else 1 - out_cubic(pr - rng.uniform(0, 0.35))
            disc(c, cx + rr * math.cos(ang) + 4 * math.sin(t * 1.2 + j), 600 + rr * math.sin(ang) + 4 * math.cos(t * 1.1 + j),
                 7, PAPER, ca * (0.25 + 0.7 * fade_out))
        text(c, lab, cx, 820, 36, PAPER, alpha=ca)
        rrect(c, cx - 110, 870, 220, 16, 8, PAPER, ca * 0.15)
        w = 220 * lerp(1.0, final, out_cubic(pr))
        rrect(c, cx - 110, 870, max(w, 6), 16, 8, col, ca)
    # the reach of the edit
    for k, cx in enumerate((1130, 1560)):
        p = clamp((t - 1.9 - k * 0.15) / 0.9)
        if 0 < p < 1:
            ring(c, 330 + (cx - 330 - 160) * out_cubic(p), 600, 20 + 50 * p, CORAL, 4, 1 - p)
    a, dy = rise(t, 3.6, 0.9)
    text(c, "goal: lower in the target cells, steady in the rest", W / 2, 990, 38, MINT, alpha=a, bold=False, dy=dy)


def icon_splice(c, x, y, p, col, dark):
    rrect(c, x - 110, y - 10, 70, 44, 10, col, 0.85 * p)
    rrect(c, x + 40, y - 10, 70, 44, 10, col, 0.85 * p)
    c.set_source_rgba(*INK, 0.8 * p)
    c.set_line_width(4)
    c.move_to(x - 75, y - 10)
    c.curve_to(x - 45, y - 70, x + 45, y - 70, x + 75, y - 10)
    c.stroke()


def icon_methyl(c, x, y, p, col, dark, t):
    line(c, x - 120, y + 20, x + 120, y + 20, 6, INK, 0.7 * p)
    for i in range(7):
        xx = x - 100 + i * 33
        bob = 3 * math.sin(t * 2 + i)
        line(c, xx, y + 20, xx, y - 12 + bob, 2, INK, 0.4 * p)
        disc(c, xx, y - 22 + bob, 11, col, p)


def icon_neighbours(c, x, y, p, col, dark, t):
    for k, dx in enumerate((-110, 0, 110)):
        rrect(c, x + dx - 40, y - 18, 80, 56, 12, col if k == 1 else INK, (0.95 if k == 1 else 0.28) * p)
    r = (t % 2.4) / 2.4
    ring(c, x, y + 10, 20 + 55 * r, col, 3, (1 - r) * p)


def icon_bystander(c, x, y, p, col, dark, t):
    s = "ATCGGCA"
    for i, ch in enumerate(s):
        xx = x - 105 + i * 35
        hit = i in (3, 4)
        text(c, ch, xx, y + 16, 40, CORAL if hit else INK, alpha=p * (1 if hit else 0.55), family=MONO)
    c.set_source_rgba(*INK, 0.6 * p)
    c.set_line_width(3)
    c.rectangle(x - 105 + 2 * 35 - 14, y - 22, 4 * 35, 70)
    c.stroke()


def icon_tissues(c, x, y, p, col, dark, t):
    cols = [(0.62, 0.55, 0.80), (0.40, 0.65, 0.80), (0.92, 0.50, 0.50), (0.55, 0.78, 0.55), GOLD]
    for i, cc in enumerate(cols):
        disc(c, x - 104 + i * 52, y + 16 + 4 * math.sin(t * 1.6 + i), 20, cc, p)


def icon_cells(c, x, y, p, col, dark, t):
    for dx in (-60, 60):
        ring(c, x + dx, y + 14, 42, col, 5, p)
        disc(c, x + dx - 8, y + 8, 14, col, 0.4 * p)


def scene_stress(c, t, dark):
    header(c, t, 3, "Stress-test every edit", dark, "Six checks, with limits fixed before any result.")
    cards = [("Splicing", icon_splice, True), ("Methylation", icon_methyl, False), ("Neighbouring genes", icon_neighbours, True),
             ("Bystander edits", icon_bystander, True), ("Other tissues", icon_tissues, True), ("Healthy cells", icon_cells, True)]
    pw, ph, gx, gy = 500, 290, 70, 40
    x0 = (W - (3 * pw + 2 * gx)) / 2
    for i, (name, fn, ok) in enumerate(cards):
        col_i, row = i % 3, i // 3
        x, y = x0 + col_i * (pw + gx), 300 + row * (ph + gy)
        a, dy = rise(t, 0.6 + i * 0.35, 0.9, 40)
        if a <= 0:
            continue
        bad = (not ok) and t > 3.4 + i * 0.1
        rrect(c, x, y + dy, pw, ph, 30, CARD, a)
        if bad:
            rrect(c, x, y + dy, pw, ph, 30, CORAL, a * 0.9, stroke=5)
        c.save()
        c.translate(0, dy)
        args = (c, x + pw / 2, y + 110, a, (CORAL if bad else TEAL), dark)
        if fn in (icon_methyl, icon_neighbours, icon_bystander, icon_tissues, icon_cells):
            fn(*args, t)
        else:
            fn(*args)
        c.restore()
        text(c, name, x + 34, y + ph - 52 + dy, 38, INK, align="left", alpha=a, family=SERIF, bold=False)
        res = clamp((t - (2.6 + i * 0.25)) / 0.5)
        if res > 0:
            if ok:
                tick(c, x + pw - 60, y + ph - 52 + dy, out_cubic(res), TEAL)
            else:
                cross(c, x + pw - 60, y + ph - 52 + dy, out_cubic((t - 3.4) / 0.5), CORAL)
    a, dy = rise(t, 4.0, 0.8)
    text(c, "most edits trip at least one wire", W / 2, 1010, 36, MUTE_L, alpha=a, bold=False, dy=dy)


def scene_judge(c, t, dark):
    if t < 3.5:
        a = 1 - fade(t, 3.0, 0.6)
        c.push_group()
        header(c, t, 4, "Judge on patients it never saw", dark, "One group of patients proposes. A different group decides.")
        rngA, rngB = np.random.default_rng(1), np.random.default_rng(2)
        for grp, (cx, col, lab, rng) in enumerate(((560, MINT, "Group A  ·  proposes", rngA), (1360, GOLD, "Group B  ·  judges", rngB))):
            ga = fade(t, 0.9 + grp * 0.4)
            for j in range(18):
                ang, rr = rng.uniform(0, 2 * math.pi), 170 * math.sqrt(rng.uniform(0.05, 1))
                disc(c, cx + rr * math.cos(ang) + 6 * math.sin(t * 0.9 + j), 640 + rr * math.sin(ang) + 6 * math.cos(t * 0.8 + j), 17, col, ga * 0.9)
            ring(c, cx, 640, 210, col, 3, ga * 0.4)
            text(c, lab, cx, 920, 38, PAPER, alpha=ga)
        wa = fade(t, 1.8)
        line(c, 960, 430, 960, 850, 5, CORAL, wa)
        text(c, "never mixed", 960, 400, 32, CORAL, alpha=wa, bold=False)
        c.pop_group_to_source()
        c.paint_with_alpha(a)
    else:
        tt = t - 3.5
        header(c, tt, 5, "…and beat chance", dark, "40 random edits set the bar every candidate has to clear.")
        ax0, ax1, ay = 200, 1720, 760
        line(c, ax0, ay, ax1, ay, 4, MUTE_D, fade(tt, 0.4))
        xmax = max(DATA["cands"]) * 1.1
        def px(v):
            return ax0 + 40 + math.sqrt(min(v, xmax) / xmax) * (ax1 - ax0 - 80)
        for i, v in enumerate(DATA["null"]):
            ti = 1.0 + i * 0.03
            a, dy = rise(tt, ti, 0.6, 30)
            disc(c, px(v), ay - 30 - 26 * (i % 5) + dy * 0, 13, MINT, a * 0.75)
        cx = px(DATA["cutoff"])
        a = fade(tt, 2.2)
        line(c, cx, ay + 30, cx, ay - 260, 4, GOLD, a)
        text(c, "best random edit", cx + 16, ay - 280, 30, GOLD, align="left", alpha=a, bold=False)
        for i, v in enumerate(DATA["cands"]):
            ti = 2.8 + i * 0.12
            a, dy = rise(tt, ti, 0.7, 30)
            disc(c, px(v), ay - 40 - 38 * (i % 3), 20, CORAL, a)
        text(c, "predicted effect on the gene", ax1, ay + 70, 30, MUTE_D, align="right", alpha=fade(tt, 0.6), bold=False)


def scene_stats(c, t, dark):
    a, dy = rise(t, 0.2, 1.0)
    text(c, "What the first run found", 150, 150, 70, INK, align="left", alpha=a, family=SERIF, bold=False, dy=dy)
    tiles = [(150, 420, 62, "patients", 0.7, INK), (1030, 420, 12738, "T cells", 1.2, INK),
             (150, 740, DATA["n_edits"], "edits tested", 1.7, INK), (1030, 740, 0, "passed every check", 2.2, CORAL)]
    for x, y, v, lab, t0, col in tiles:
        a, dy = rise(t, t0, 1.0, 30)
        text(c, f"{v * out_cubic((t - t0) / 1.2):,.0f}", x, y, 150, col, align="left", alpha=a, family=SERIF, bold=False, dy=dy, fit=False)
        text(c, lab, x + 4, y + 100, 38, MUTE_L, align="left", alpha=a, bold=False, dy=dy, fit=False)
    a, dy = rise(t, 3.3, 1.0)
    text(c, f"{DATA['n_fail_meth']} of {DATA['n_edits']} were stopped by methylation.", 150, 960, 46, INK,
         align="left", alpha=a, family=SERIF, bold=False, dy=dy)


def scene_byo(c, t, dark):
    a, dy = rise(t, 0.2, 1.0)
    text(c, "Use it on your own data", 150, 150, 78, PAPER, align="left", alpha=a, family=SERIF, bold=False, dy=dy)
    a = fade(t, 0.5)
    rrect(c, 150, 260, 1000, 600, 28, (0.09, 0.19, 0.17), a)
    for k, col in enumerate((CORAL, GOLD, MINT)):
        disc(c, 200 + k * 36, 304, 9, col, a)
    lines = ["$ pip install gene-edit-check", "$ gene-edit-check doctor", "$ gene-edit-check run-h5ad \\",
             "    --h5ad my_cells.h5ad \\", "    --gene-symbol PDCD1 \\", "    --scan-region 14300:15300"]
    total = sum(len(s) for s in lines)
    left = int(total * clamp((t - 0.9) / 2.6))
    y = 380
    for s in lines:
        n = min(len(s), left)
        left -= n
        if n <= 0:
            break
        indent = len(s) - len(s.lstrip())
        text(c, s[:n].lstrip() + ("_" if left <= 0 and int(t * 2) % 2 == 0 else ""), 200 + indent * 21, y, 38,
             MINT if s.startswith("$") else PAPER, align="left", family=MONO, bold=False, fit=False)
        y += 78
    chips = [("any gene", 1.2), ("your own single-cell data", 1.7), ("AI agents through MCP", 2.2), ("macOS · Linux · Windows", 2.7)]
    for i, (s, t0) in enumerate(chips):
        a, dy = rise(t, t0, 0.9, 30)
        rrect(c, 1260, 280 + i * 140, 520, 96, 48, PAPER, a * 0.1)
        rrect(c, 1260, 280 + i * 140, 520, 96, 48, MINT, a * 0.7, stroke=2)
        text(c, s, 1520, 328 + i * 140, 32, PAPER, alpha=a, dy=dy)
    text(c, "Any gene. Any cells. One command.", W / 2, 960, 52, MINT, alpha=fade(t, 3.4), family=SERIF, bold=False)


def scene_end(c, t, dark):
    helix(c, 470, 540, 920, 160, 0.8 + t * 0.5, grow=out_cubic(t / 1.6), dot=14, n=24)
    a, dy = rise(t, 0.6, 1.2)
    w = text(c, "gene-edit-", 900, 400, 110, PAPER, align="left", alpha=a, family=SERIF, bold=False, dy=dy)
    text(c, "check", 900 + w, 400, 110, MINT, align="left", alpha=a, family=SERIF, bold=False, dy=dy)
    a, dy = rise(t, 1.4, 1.0)
    text(c, "Check a gene edit on a computer", 904, 520, 42, PAPER, align="left", alpha=a, bold=False, dy=dy)
    text(c, "before anyone spends weeks in a lab.", 904, 582, 42, PAPER, align="left", alpha=a, bold=False, dy=dy)
    s = "pip install gene-edit-check"
    n = int(len(s) * clamp((t - 2.2) / 1.2))
    text(c, s[:n] + ("_" if int(t * 2) % 2 == 0 else ""), 904, 730, 48, GOLD, align="left", family=MONO, bold=False)
    text(c, "github.com/aravmendiratta/gene-edit-check", 904, 820, 34, MUTE_D, align="left", alpha=fade(t, 3.4), bold=False)


SCENE_FN = {"hook": scene_hook, "find": scene_find, "predict": scene_predict, "stress": scene_stress,
            "judge": scene_judge, "stats": scene_stats, "byo": scene_byo, "end": scene_end}


def draw_scene(c, idx, t):
    name, _, dark = SCENES[idx]
    c.set_source_rgb(*(FOREST if dark else PAPER))
    c.paint()
    c.save()
    z = 1.0 + 0.02 * clamp(t / (SCENES[idx][1] * BAR))  # a very slow push-in
    c.translate(W / 2, H / 2)
    c.scale(z, z)
    c.translate(-W / 2, -H / 2)
    SCENE_FN[name](c, t, dark)
    c.restore()


def render(t, surf):
    c = cairo.Context(surf)
    c.scale(OW / W, OH / H)
    idx = max(i for i in range(len(SCENES)) if STARTS[i] <= t + 1e-9)
    local = t - STARTS[idx]
    if idx > 0 and local < 0.7:
        draw_scene(c, idx - 1, SCENES[idx - 1][1] * BAR + local)
        c.push_group()
        draw_scene(c, idx, local)
        c.pop_group_to_source()
        c.paint_with_alpha(in_out(local / 0.7))
    else:
        draw_scene(c, idx, local)


# ---- music -----------------------------------------------------------------------------------
SR = 44100


def make_music() -> np.ndarray:
    """Soft and unhurried: warm pads, celesta-like bells, a low pulse that arrives late. 100 BPM, A minor."""
    from scipy.signal import butter, fftconvolve, lfilter

    n = int(TOTAL * SR)
    dry = np.zeros((n, 2))
    wet = np.zeros((n, 2))
    rng = np.random.default_rng(9)

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

    def bell(f, d=2.0):
        t = tt(d)
        s = np.sin(2 * np.pi * f * t) + 0.35 * np.sin(2 * np.pi * 2 * f * t) * np.exp(-t * 4) + 0.15 * np.sin(2 * np.pi * 3.01 * f * t) * np.exp(-t * 7)
        return s * np.exp(-t * 2.4) * np.minimum(1, t * 200) * 0.30

    def thump(d=0.4):
        t = tt(d)
        return np.sin(2 * np.pi * np.cumsum(52 + 60 * np.exp(-t * 30)) / SR) * np.exp(-t * 9) * 0.55

    chords = [(110.0, [220.0, 261.63, 329.63]), (87.31, [174.61, 220.0, 261.63]),
              (130.81, [196.0, 261.63, 329.63]), (98.0, [196.0, 246.94, 293.66])]
    bars = int(round(TOTAL / BAR))
    for bar in range(bars):
        root, triad = chords[bar % 4]
        t_bar = bar * BAR
        d = BAR + 1.4
        env = np.minimum(1, tt(d) / 0.9) * np.minimum(1, (d - tt(d)) / 1.2)
        sl = sum(saw(f, d, de) for f in triad + [triad[0] * 2] for de in (-0.003, 0.002)) / 8
        sr_ = sum(saw(f, d, de) for f in triad + [triad[0] * 2] for de in (0.003, -0.002)) / 8
        put(wet, t_bar, filt(sl, "low", 1100) * env * 0.55, 1.0, 0.0)
        put(wet, t_bar, filt(sr_, "low", 1100) * env * 0.55, 0.0, 1.0)
        if bar >= 3:
            put(dry, t_bar, np.sin(2 * np.pi * root * tt(BAR * 1.1)) * np.exp(-tt(BAR * 1.1) * 1.1) * 0.32)
        # slow bells: a few notes per bar, drawn from the chord, answering each other left to right
        pattern = [0, 2.0, 3.0] if bar % 2 == 0 else [0.5, 1.5, 3.0]
        for k, b in enumerate(pattern):
            if bar < 1:
                break
            f = triad[(bar + k) % 3] * (2 if (k % 2 == 0) else 4)
            tb = t_bar + b * BEAT
            pan = 0.35 + 0.3 * (k % 2)
            put(wet, tb, bell(f), 1 - pan, pan)
            put(dry, tb, bell(f) * 0.25, 1 - pan, pan)
            put(wet, tb + BEAT * 0.75, bell(f) * 0.35, pan, 1 - pan)
        if bar >= 4:
            put(dry, t_bar, thump(), 0.8)
            put(dry, t_bar + 2 * BEAT, thump(), 0.55)
        if bar >= 7:
            for e in range(8):
                put(dry, t_bar + e * BEAT / 2, filt(rng.normal(0, 1, int(0.05 * SR)), "high", 6500) * np.exp(-tt(0.05) * 70) * 0.05, 0.8, 0.9)
    # a gentle swell into every scene change, and one soft low note on the cut
    for s in range(1, len(SCENES)):
        t0 = float(STARTS[s])
        d = 1.6
        sw = filt(rng.normal(0, 1, int(d * SR)), "low", 2200) * (tt(d) / d) ** 2 * 0.10
        put(wet, t0 - d, sw)
        put(dry, t0, np.sin(2 * np.pi * 55 * tt(1.8)) * np.exp(-tt(1.8) * 2.2) * 0.30)
    # reverb: a long, dark tail
    ir_t = tt(2.6)
    ir = filt(rng.normal(0, 1, len(ir_t)), "low", 3500) * np.exp(-ir_t * 2.4)
    ir /= np.abs(ir).sum() / 14
    rev = np.stack([fftconvolve(wet[:, 0], ir)[:n], fftconvolve(wet[:, 1], np.roll(ir, 211))[:n]], axis=1)
    mix = dry + rev * 0.9 + wet * 0.35
    fade_in = np.minimum(1, np.arange(n) / SR / 1.0)
    fade_out = np.minimum(1, (n - np.arange(n)) / SR / 2.4)
    mix *= (fade_in * fade_out)[:, None]
    mix = np.tanh(mix * 0.9)
    mix /= max(1e-6, np.abs(mix).max()) / 0.70
    return (mix * 32767).astype(np.int16)


# ---- encode ----------------------------------------------------------------------------------
def main(out: str = "assets/launch.mp4") -> None:
    import wave

    import imageio_ffmpeg

    DATA.update(load_data())
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    wav = Path(out).with_suffix(".wav")
    audio = make_music()
    with wave.open(str(wav), "wb") as w:
        w.setnchannels(2)
        w.setsampwidth(2)
        w.setframerate(SR)
        w.writeframes(audio.tobytes())
    cmd = [imageio_ffmpeg.get_ffmpeg_exe(), "-y", "-loglevel", "error",
           "-f", "rawvideo", "-pix_fmt", "bgra", "-s", f"{OW}x{OH}", "-r", str(FPS), "-i", "-", "-i", str(wav),
           "-c:v", "libx264", "-profile:v", "main", "-level", "4.0", "-preset", "slow", "-crf", "22",
           "-maxrate", "3M", "-bufsize", "6M", "-bf", "0", "-g", str(FPS), "-keyint_min", str(FPS), "-sc_threshold", "0",
           "-pix_fmt", "yuv420p", "-movflags", "+faststart",
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
