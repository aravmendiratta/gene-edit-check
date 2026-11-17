"""A fast 30-second motion-graphics launch video for the tool, drawn entirely in code.

    pip install pycairo imageio imageio-ffmpeg numpy scipy
    PYTHONPATH=src python scripts/make_launch_video.py        # writes assets/launch.mp4

No screenshots and no stock footage. The numbers on screen (the candidate bars, the random-edit null, the
patient and cell counts) are read from the saved results in `outputs/`. The music is synthesised here
(128 BPM, A minor), so there is nothing to license. Video and sound share one ffmpeg pass with the index at the
front of the file, which is what makes it play without stalling.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import cairo
import numpy as np

W, H, FPS = 1920, 1080, 30
BPM = 128
BEAT = 60 / BPM
BAR = BEAT * 4
SCENE_LEN = 2 * BAR  # 3.75 s: every scene is two bars, so cuts land on the music
OUT = Path("outputs")

INK = (0.043, 0.043, 0.071)
YELLOW = (1.0, 0.84, 0.04)
PINK = (1.0, 0.24, 0.51)
CYAN = (0.18, 0.89, 0.90)
BLUE = (0.17, 0.29, 1.0)
WHITE = (1, 1, 1)
GREY = (0.45, 0.46, 0.52)
PANEL = (0.10, 0.10, 0.17)

SCENES = [  # name, background
    ("hook", INK), ("propose", INK), ("predict", BLUE), ("stress", INK),
    ("honest", YELLOW), ("stats", PINK), ("byo", INK), ("end", INK),
]
TOTAL = SCENE_LEN * len(SCENES)
DATA: dict = {}


# ---- easing and helpers ----------------------------------------------------------------------
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


FONT = "Arial Black"


def set_font(c: cairo.Context, size: float, bold: bool = True, family: str | None = None) -> None:
    c.select_font_face(family or FONT, cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD if bold else cairo.FONT_WEIGHT_NORMAL)
    c.set_font_size(size)


def text(c: cairo.Context, s: str, x: float, y: float, size: float, color: tuple, *, align: str = "center",
         scale: float = 1.0, scale_y: float | None = None, alpha: float = 1.0, rot: float = 0.0,
         family: str | None = None, bold: bool = True, fit: bool = True) -> float:
    """Draw text centred on y. Long lines shrink to stay on screen. Returns the drawn width."""
    if alpha <= 0.001 or scale <= 0.001:
        return 0.0
    c.save()
    set_font(c, size, bold, family)
    ext = c.text_extents(s)
    if fit and ext.width > 0:
        avail = {"left": W - x - 70, "right": x - 70}.get(align, 2 * min(x, W - x) - 100)
        if ext.width * scale > avail and scale <= 1.05:
            scale *= avail / (ext.width * scale)
    c.translate(x, y)
    c.rotate(rot)
    c.scale(scale, scale if scale_y is None else scale_y)
    ox = {"center": -ext.width / 2 - ext.x_bearing, "left": -ext.x_bearing, "right": -ext.width - ext.x_bearing}[align]
    c.move_to(ox, size * 0.36)
    c.set_source_rgba(*color, alpha)
    c.show_text(s)
    c.restore()
    return ext.width * scale


def slam(t: float, t0: float, dur: float = 0.26) -> tuple[float, float]:
    if t < t0:
        return 0.0, 0.0
    p = (t - t0) / dur
    return lerp(1.8, 1.0, out_back(p, 1.2)), clamp(p * 4)


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


def rrect(c: cairo.Context, x: float, y: float, w: float, h: float, r: float, color: tuple, alpha: float = 1.0) -> None:
    c.set_source_rgba(*color, alpha)
    c.new_sub_path()
    c.arc(x + w - r, y + r, r, -math.pi / 2, 0)
    c.arc(x + w - r, y + h - r, r, 0, math.pi / 2)
    c.arc(x + r, y + h - r, r, math.pi / 2, math.pi)
    c.arc(x + r, y + r, r, math.pi, 3 * math.pi / 2)
    c.close_path()
    c.fill()


def tick(c: cairo.Context, x: float, y: float, p: float, color: tuple, w: float = 11) -> None:
    c.set_source_rgba(*color, 1)
    c.set_line_width(w * p)
    c.set_line_cap(cairo.LINE_CAP_ROUND)
    c.move_to(x - 24 * p, y + 2 * p)
    c.line_to(x - 6 * p, y + 22 * p)
    c.line_to(x + 28 * p, y - 22 * p)
    c.stroke()


def cross(c: cairo.Context, x: float, y: float, p: float, color: tuple, w: float = 12) -> None:
    c.set_source_rgba(*color, 1)
    c.set_line_width(w * p)
    c.set_line_cap(cairo.LINE_CAP_ROUND)
    c.move_to(x - 22 * p, y - 22 * p)
    c.line_to(x + 22 * p, y + 22 * p)
    c.move_to(x + 22 * p, y - 22 * p)
    c.line_to(x - 22 * p, y + 22 * p)
    c.stroke()


BASE_COLORS = {"A": CYAN, "T": YELLOW, "C": PINK, "G": WHITE}
PAIR = {"A": "T", "T": "A", "C": "G", "G": "C"}


def helix(c: cairo.Context, cx: float, cy: float, height: float, amp: float, rot: float, *, alpha: float = 1.0,
          n: int = 26, turns: float = 2.4, grow: float = 1.0, seq: str = "ATCGGCTAACGTTGCAGCTA", dot: float = 11.0) -> None:
    items = []
    for i in range(n):
        u = i / (n - 1)
        if u > grow:
            break
        y = cy - height / 2 + height * u
        ph = u * turns * 2 * math.pi + rot
        xa, xb = cx + amp * math.sin(ph), cx - amp * math.sin(ph)
        items.append((i, y, xa, xb, math.cos(ph), -math.cos(ph), seq[i % len(seq)]))
    for i, y, xa, xb, da, db, b in items:
        line(c, xa, y, (xa + xb) / 2, y, 4, BASE_COLORS[b], alpha * 0.9)
        line(c, (xa + xb) / 2, y, xb, y, 4, BASE_COLORS[PAIR[b]], alpha * 0.9)
    for k in range(len(items) - 1):
        _, y0, xa0, xb0, da0, db0, _ = items[k]
        _, y1, xa1, xb1, _, _, _ = items[k + 1]
        line(c, xa0, y0, xa1, y1, 5, WHITE, alpha * (0.22 + 0.3 * max(da0, 0)))
        line(c, xb0, y0, xb1, y1, 5, WHITE, alpha * (0.22 + 0.3 * max(db0, 0)))
    dots = []
    for i, y, xa, xb, da, db, b in items:
        dots += [(da, xa, y, BASE_COLORS[b]), (db, xb, y, BASE_COLORS[PAIR[b]])]
    for d, x, y, col in sorted(dots, key=lambda q: q[0]):
        disc(c, x, y, dot * (0.75 + 0.35 * d), col, alpha * (0.55 + 0.45 * (d + 1) / 2))


def counter(t: float, t0: float, dur: float, to: float) -> float:
    return to * out_cubic((t - t0) / dur)


# ---- data ------------------------------------------------------------------------------------
def load_data() -> dict:
    rep = json.loads((OUT / "report_crc-pd1-exhaustion.json").read_text(encoding="utf-8"))
    its = rep["iterations"]
    cands = []
    for it in its:
        t1 = it["report"]["tier1_selectivity"]
        h = t1.get("healthy_control_deltas") or {}
        hv = next(iter(h.values())) if isinstance(h, dict) and h else 0.0
        if isinstance(hv, dict):
            hv = next(iter(hv.values()), 0.0)
        cands.append({"pos": it["report"]["candidate"]["target_edit"]["position"], "lfc": t1["target_delta_sketch_b"],
                      "healthy": float(hv) if isinstance(hv, (int, float)) else 0.0})
    null = rep["efficacy_null"]
    return {"cands": cands, "null": null["abs_lfc_values"], "cutoff": null["abs_lfc_max"], "n_edits": len(its),
            "n_fail_meth": sum(1 for it in its if not it["report"]["tier2_hazards"]["methylation_hazard_passed"])}


# ---- scenes (t is seconds since the scene started) --------------------------------------------
def scene_hook(c, t):
    """Chaos first, then the name lands."""
    fr = int(t * FPS)
    rng = np.random.default_rng(fr)
    chaos = t < 2.95
    if chaos:
        bg = [INK, INK, INK, YELLOW, PINK, WHITE][int(rng.integers(0, 6)) if (fr % 3 == 0) else 0]
        if t < 0.1:
            bg = WHITE
        c.set_source_rgb(*bg)
        c.paint()
        c.translate(rng.normal(0, 24 * (1 - t / 3.0)), rng.normal(0, 24 * (1 - t / 3.0)))
    z = 1 + 0.06 * math.sin(t * 18) * (1 - t / 3.0)
    c.translate(W / 2, H / 2)
    c.scale(z, z)
    c.translate(-W / 2, -H / 2)
    if chaos:
        helix(c, W / 2, H / 2, 1180, 260, t * (2 + t * 5), alpha=0.9, n=30, dot=19)
        g = np.random.default_rng(7)
        for k in range(150):
            s0 = g.uniform(0.05, 2.5)
            ang = g.uniform(0, 2 * math.pi)
            sp = g.uniform(350, 1500)
            size = 60 + 120 * g.uniform()
            spin = g.uniform(-4, 4)
            if t < s0:
                continue
            age = t - s0
            r = sp * age ** 1.25
            x, y = W / 2 + math.cos(ang) * r, H / 2 + math.sin(ang) * r
            if -100 < x < W + 100 and -100 < y < H + 100:
                ch = "ATCG"[k % 4]
                text(c, ch, x, y, size, BASE_COLORS[ch], rot=ang + age * spin, alpha=clamp(1 - age / 1.6), fit=False)
        words = ["GENES", "EDITS", "CELLS", "RISK", "WEEKS", "LAB", "???", "DATA", "ERROR"]
        if 0.3 < t < 2.4:
            k = int((t - 0.3) / 0.16)
            r2 = np.random.default_rng(100 + k)
            col = [WHITE, YELLOW, PINK, CYAN][k % 4]
            text(c, words[k % len(words)], r2.uniform(500, 1400), r2.uniform(330, 760), r2.uniform(260, 470), col,
                 rot=r2.uniform(-0.18, 0.18), fit=False)
        if 2.45 <= t < 2.95:
            sc, al = slam(t, 2.45, 0.15)
            text(c, "CHECK IT", W / 2, 400, 330, WHITE, scale=sc, alpha=al)
            sc, al = slam(t, 2.7, 0.15)
            text(c, "FIRST.", W / 2, 720, 330, PINK, scale=sc, alpha=al)
    else:
        c.identity_matrix()
        c.set_source_rgb(*INK)
        c.paint()
        helix(c, W / 2, H / 2, 1100, 240, 8 + t * 2, alpha=0.28, n=30, dot=17)
        sc, al = slam(t, 2.97, 0.3)
        c.save()
        c.translate(W / 2, 520)
        c.scale(max(sc, 0.001), max(sc, 0.001))
        set_font(c, 170)
        e1, e2 = c.text_extents("gene-edit-"), c.text_extents("check")
        x0 = -(e1.x_advance + e2.x_advance) / 2
        c.move_to(x0, 170 * 0.36)
        c.set_source_rgba(*WHITE, al)
        c.show_text("gene-edit-")
        c.set_source_rgba(*CYAN, al)
        c.show_text("check")
        c.restore()
        flash = clamp(1 - (t - 2.95) / 0.25)
        if flash > 0:
            c.set_source_rgba(1, 1, 1, flash * 0.9)
            c.paint()


def scene_propose(c, t):
    sc, al = slam(t, 0.05)
    text(c, "PICK A GENE.", 120, 130, 110, WHITE, align="left", scale=sc, alpha=al)
    sc, al = slam(t, 0.45)
    text(c, "IT FINDS EVERY EDIT.", 120, 270, 110, YELLOW, align="left", scale=sc, alpha=al)
    seq = "ACGTTGGCATCGGATCCGTAAGGCTTCAGGACCTGGATGCAAGG"
    x0, step, y = 130, 38, 720
    idx = (t - 0.5) / 2.6 * len(seq)
    found = 0
    for i, ch in enumerate(seq):
        pop = out_back((t - 0.1 - i * 0.012) / 0.25, 2.0)
        text(c, ch, x0 + i * step, y, 54 * pop, BASE_COLORS[ch] if i < idx else GREY, family="Consolas", fit=False)
    for i in range(1, len(seq) - 1):
        if seq[i] == "G" and seq[i + 1] == "G" and idx > i + 1:
            found += 1
            p = out_back((idx - i - 1) / 1.2, 1.6)
            col = PINK if found % 2 else CYAN
            c.set_source_rgba(*col, 0.9)
            c.set_line_width(5)
            c.rectangle(x0 + (i - 1) * step - 17, y - 42, step * 3, 84 * p)
            c.stroke()
            text(c, "ABE8e" if found % 2 else "CBE", x0 + i * step, y - 80, 30, col, alpha=clamp(p), family="Arial", fit=False)
    scan = x0 + idx * step
    if 0 <= idx <= len(seq):
        c.set_source_rgba(*WHITE, 0.55)
        c.rectangle(scan - 6, y - 60, 12, 120)
        c.fill()
    v = counter(t, 0.6, 2.4, 60)
    text(c, f"{v:.0f}", 1560, 470, 240, PINK)
    text(c, "CANDIDATES RANKED", 1560, 590, 46, WHITE)
    sc, al = slam(t, 2.6)
    text(c, "BY WHERE THE MODEL IS MOST SENSITIVE", W / 2, 940, 56, WHITE, scale=sc, alpha=al)


def scene_predict(c, t):
    sc, al = slam(t, 0.05)
    text(c, "THEN IT PREDICTS EACH ONE.", 120, 120, 100, WHITE, align="left", scale=sc, alpha=al)
    cands = DATA["cands"][:4]
    for i, cd in enumerate(cands):
        y = 330 + i * 160
        t_in = 0.5 + i * 0.3
        a = clamp((t - t_in) / 0.25)
        text(c, f"chr2:{cd['pos']:,}", 120, y + 20, 40, WHITE, align="left", family="Consolas", alpha=a, fit=False)
        full = 760
        for k, (val, col, dy, h) in enumerate(((cd["lfc"], PINK, 0, 54), (cd["healthy"], CYAN, 62, 34))):
            prog = out_cubic((t - t_in - 0.2 - k * 0.1) / 0.8)
            ln = full * (1 - prog) + full * (2 ** val) * prog
            rrect(c, 620, y - 14 + dy, max(ln, 4), h, 10, col, a)
            if prog > 0.9:
                pct = (2 ** val - 1) * 100
                text(c, f"{pct:+.0f}%", 620 + ln + 20, y - 14 + dy + h / 2, 36, WHITE, align="left",
                     alpha=clamp((prog - 0.9) * 10), family="Arial", fit=False)
    rrect(c, 1500, 330, 30, 30, 8, PINK)
    text(c, "TARGET CELLS", 1545, 345, 34, WHITE, align="left", family="Arial", fit=False)
    rrect(c, 1500, 390, 30, 30, 8, CYAN)
    text(c, "HEALTHY CELLS", 1545, 405, 34, WHITE, align="left", family="Arial", fit=False)
    sc, al = slam(t, 2.5)
    text(c, "ALPHAGENOME 2, ONE CELL TYPE AT A TIME.", W / 2, 960, 52, WHITE, scale=sc, alpha=al)


def scene_stress(c, t):
    sc, al = slam(t, 0.05)
    text(c, "THEN IT TRIES TO BREAK THEM.", 120, 120, 100, WHITE, align="left", scale=sc, alpha=al)
    names = ["SPLICING", "METHYLATION", "NEIGHBOUR GENES", "BYSTANDER EDITS", "OTHER TISSUES", "HEALTHY CELLS"]
    pw, ph, gx, gy = 540, 220, 60, 50
    x0 = (W - (3 * pw + 2 * gx)) / 2
    for i, n in enumerate(names):
        col, row = i % 3, i // 3
        x, y = x0 + col * (pw + gx), 260 + row * (ph + gy)
        t_in = 0.35 + i * 0.12
        p = out_back((t - t_in) / 0.35, 1.6)
        if p <= 0:
            continue
        yy = y + (1 - clamp(p)) * 300
        is_bad = n == "METHYLATION"
        res = clamp((t - (1.8 + i * 0.16)) / 0.25)
        showing_bad = is_bad and res > 0
        rrect(c, x, yy, pw, ph, 28, YELLOW if showing_bad else PANEL, clamp(p))
        text(c, n, x + 34, yy + ph / 2, 50, INK if showing_bad else WHITE, align="left", alpha=clamp(p), fit=False)
        if res > 0:
            if is_bad:
                cross(c, x + pw - 70, yy + ph / 2, out_back(res, 2.2), INK)
            else:
                tick(c, x + pw - 70, yy + ph / 2, out_back(res, 2.2), CYAN)
    sc, al = slam(t, 2.9)
    text(c, "SIX CHECKS. FIXED LIMITS. NO NEGOTIATING.", W / 2, 960, 62, YELLOW, scale=sc, alpha=al)


def scene_honest(c, t):
    sc, al = slam(t, 0.05)
    text(c, "NO LUCKY WINS.", 120, 130, 130, INK, align="left", scale=sc, alpha=al)
    sc, al = slam(t, 0.45)
    text(c, "40 RANDOM EDITS SET THE BAR.", 120, 260, 80, PINK, align="left", scale=sc, alpha=al)
    ax0, ax1, ay = 200, 1720, 700
    line(c, ax0, ay, ax1, ay, 6, INK)
    xmax = 0.4
    for i, v in enumerate(DATA["null"]):
        t_in = 0.8 + i * 0.025
        r = 15 * out_back((t - t_in) / 0.25, 2.0)
        disc(c, ax0 + 30 + v / xmax * (ax1 - ax0 - 60), ay - 24 - 26 * (i % 5) * out_cubic((t - t_in) / 0.4), r, INK, 0.85)
    cx = ax0 + 30 + DATA["cutoff"] / xmax * (ax1 - ax0 - 60)
    a = clamp((t - 1.9) / 0.3)
    line(c, cx, ay + 30, cx, ay - 250, 5, PINK, a)
    text(c, "BEST RANDOM EDIT", cx + 14, ay - 275, 32, PINK, align="left", alpha=a, family="Arial", fit=False)
    for i, cd in enumerate(DATA["cands"]):
        t_in = 2.1 + i * 0.12
        r = 26 * out_back((t - t_in) / 0.3, 2.2)
        px = ax0 + 30 + min(abs(cd["lfc"]), xmax) / xmax * (ax1 - ax0 - 60)
        disc(c, px, ay - 36 - 40 * (i % 3), r, PINK)
    text(c, "predicted effect on the gene →", ax1, ay + 60, 32, INK, align="right", family="Arial", bold=False, fit=False, alpha=0.8)
    sc, al = slam(t, 2.8)
    text(c, "ONE GROUP OF PATIENTS PROPOSES. A NEW GROUP JUDGES.", W / 2, 940, 50, INK, scale=sc, alpha=al)


def scene_stats(c, t):
    tiles = [(480, 330, 62, "PATIENTS", 0.1), (1440, 330, 12738, "T CELLS", 0.6),
             (480, 700, DATA["n_edits"], "EDITS TESTED", 1.1), (1440, 700, 0, "PASSED EVERY CHECK", 1.6)]
    for x, y, v, lab, t0 in tiles:
        sc, al = slam(t, t0, 0.22)
        shown = v * out_cubic((t - t0) / 0.6)
        text(c, f"{shown:,.0f}", x, y, 260, WHITE if v else INK, scale=sc, alpha=al)
        text(c, lab, x, y + 185, 52, INK, scale=sc, alpha=al)
    sc, al = slam(t, 2.6)
    text(c, f"{DATA['n_fail_meth']} OF {DATA['n_edits']} FAIL ON METHYLATION.", W / 2, 970, 80, WHITE, scale=sc, alpha=al)


def scene_byo(c, t):
    sc, al = slam(t, 0.05)
    text(c, "USE IT ON YOUR OWN DATA.", 120, 120, 100, WHITE, align="left", scale=sc, alpha=al)
    rrect(c, 120, 230, 1020, 640, 26, PANEL)
    for k, col in enumerate((PINK, YELLOW, CYAN)):
        disc(c, 170 + k * 40, 276, 11, col)
    lines = ["$ pip install gene-edit-check", "$ gene-edit-check doctor", "$ gene-edit-check optimize",
             "$ gene-edit-check run-h5ad \\", "    --h5ad my_cells.h5ad \\", "    --gene-symbol PDCD1 \\", "    --scan-region 14300:15300"]
    total = sum(len(s) for s in lines)
    left = int(total * clamp((t - 0.4) / 2.2))
    y = 340
    for s in lines:
        n = min(len(s), left)
        left -= n
        if n <= 0:
            break
        text(c, s[:n] + ("_" if left <= 0 and int(t * 4) % 2 == 0 else ""), 170, y, 42, CYAN if s.startswith("$") else WHITE,
             align="left", family="Consolas", fit=False)
        y += 70
    chips = [("ANY GENE", PINK, 0.8), ("YOUR OWN .H5AD CELLS", YELLOW, 1.1), ("AI AGENTS VIA MCP", CYAN, 1.4),
             ("MAC  ·  LINUX  ·  WINDOWS", WHITE, 1.7)]
    for i, (s, col, t0) in enumerate(chips):
        p = out_back((t - t0) / 0.35, 1.8)
        if p <= 0:
            continue
        off = (1 - clamp(p)) * 60
        rrect(c, 1230, 270 + i * 150 - off, 570, 110, 55, col, clamp(p))
        text(c, s, 1230 + 285, 325 + i * 150 - off, 36, INK, alpha=clamp(p))
    sc, al = slam(t, 2.8)
    text(c, "ANY GENE. ANY CELLS. ONE COMMAND.", W / 2, 960, 70, WHITE, scale=sc, alpha=al)


def scene_end(c, t):
    helix(c, 470, 540, 940, 165, t * 1.6 + 2, grow=out_cubic(t / 1.0), dot=15, n=24)
    sc, al = slam(t, 0.3)
    w = text(c, "gene-edit-", 880, 380, 104, WHITE, align="left", scale=sc, alpha=al)
    text(c, "check", 880 + w, 380, 104, CYAN, align="left", scale=sc, alpha=al)
    sc, al = slam(t, 0.8)
    text(c, "Check a gene edit on a computer", 880, 510, 46, WHITE, align="left", scale=sc, alpha=al, family="Arial")
    text(c, "before anyone spends weeks in a lab.", 880, 575, 46, WHITE, align="left", scale=sc, alpha=al, family="Arial")
    s = "pip install gene-edit-check"
    n = int(len(s) * clamp((t - 1.3) / 0.8))
    text(c, s[:n] + ("_" if int(t * 3) % 2 == 0 else ""), 880, 740, 54, YELLOW, align="left", family="Consolas")
    sc, al = slam(t, 2.3)
    text(c, "github.com/aravmendiratta/gene-edit-check", 880, 840, 38, GREY, align="left", scale=sc, alpha=al, family="Arial", bold=False)


SCENE_FN = {"hook": scene_hook, "propose": scene_propose, "predict": scene_predict, "stress": scene_stress,
            "honest": scene_honest, "stats": scene_stats, "byo": scene_byo, "end": scene_end}


def draw_scene(c: cairo.Context, idx: int, t: float) -> None:
    name, bg = SCENES[idx]
    c.set_source_rgb(*bg)
    c.paint()
    c.save()
    if name != "hook":
        z = 1.0 + 0.035 * clamp(t / SCENE_LEN)
        c.translate(W / 2, H / 2)
        c.scale(z, z)
        c.translate(-W / 2, -H / 2)
    SCENE_FN[name](c, t)
    c.restore()


def render(t: float, surf: cairo.ImageSurface) -> None:
    c = cairo.Context(surf)
    idx = min(int(t / SCENE_LEN + 1e-9), len(SCENES) - 1)
    local = t - idx * SCENE_LEN
    if idx > 0 and local < 0.28:
        draw_scene(c, idx - 1, SCENE_LEN + local)
        c.save()
        c.rectangle(0, 0, W * in_out(local / 0.28), H)
        c.clip()
        draw_scene(c, idx, local)
        c.restore()
    else:
        draw_scene(c, idx, local)
    if idx > 0 and local < 0.1:
        c.set_source_rgba(1, 1, 1, 0.25 * (1 - local / 0.1))
        c.paint()


def glitch_level(t: float) -> float:
    if t < 2.95:
        return 0.9
    if t < 3.3:
        return 0.7 * (1 - (t - 2.95) / 0.35)
    local = t % SCENE_LEN
    if local < 0.25 and t >= SCENE_LEN:
        return 0.45 * (1 - local / 0.25)
    return 0.0


def glitch(arr: np.ndarray, level: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = arr.copy()
    for _ in range(int(3 + 12 * level)):
        y0 = int(rng.integers(0, H - 100))
        h = int(rng.integers(6, 30 + int(110 * level)))
        out[y0 : y0 + h] = np.roll(out[y0 : y0 + h], int(rng.normal(0, 140 * level)), axis=1)
    sh = int(14 * level)
    out[:, :, 0] = np.roll(out[:, :, 0], sh, axis=1)
    out[:, :, 2] = np.roll(out[:, :, 2], -sh, axis=1)
    return out


# ---- music -----------------------------------------------------------------------------------
SR = 44100


def make_music() -> np.ndarray:
    """Stereo int16, 128 BPM, A minor: a chaotic build, then a four-on-the-floor drop with bass, pad and arp."""
    from scipy.signal import butter, lfilter

    n = int(TOTAL * SR)
    left = np.zeros(n)
    right = np.zeros(n)
    rng = np.random.default_rng(4)

    def filt(x, kind, fc, order=2):
        b, a = butter(order, fc / (SR / 2), btype=kind)
        return lfilter(b, a, x)

    def put(t0, sig, gl=1.0, gr=None):
        i = int(t0 * SR)
        if i < 0 or i >= n:
            return
        m = min(len(sig), n - i)
        left[i : i + m] += sig[:m] * gl
        right[i : i + m] += sig[:m] * (gl if gr is None else gr)

    def tt(d):
        return np.arange(int(d * SR)) / SR

    def kick(d=0.38):
        t = tt(d)
        ph = np.cumsum(46 + 150 * np.exp(-t * 32)) * 2 * np.pi / SR
        return np.sin(ph) * np.exp(-t * 8.5) + 0.35 * rng.normal(0, 1, len(t)) * np.exp(-t * 140)

    def clap(d=0.24):
        t = tt(d)
        env = np.exp(-t * 20) * (1 + 0.7 * (np.sin(t * 2 * np.pi * 95) > 0.2))
        return filt(rng.normal(0, 1, len(t)), "high", 1400) * env * 0.5

    def hat(d=0.05, open_=False):
        d = 0.2 if open_ else d
        t = tt(d)
        return filt(rng.normal(0, 1, len(t)), "high", 7000) * np.exp(-t * (16 if open_ else 80)) * 0.28

    def saw(f, d, detune=0.0):
        return 2 * ((tt(d) * f * (1 + detune)) % 1.0) - 1

    def pluck(f, d=0.2):
        t = tt(d)
        s = saw(f, d) * 0.5 + saw(f * 1.005, d) * 0.5
        return filt(s, "low", 3800) * np.exp(-t * 13) * 0.32

    def boom(d=1.2):
        t = tt(d)
        return (np.sin(2 * np.pi * np.cumsum(70 * np.exp(-t * 2.0) + 28) / SR) * np.exp(-t * 2.6)
                + filt(rng.normal(0, 1, len(t)), "high", 3000) * np.exp(-t * 5) * 0.35)

    def riser(d):
        t = tt(d)
        s = filt(rng.normal(0, 1, len(t)), "high", 1200) * (t / d) ** 2.2 * 0.3
        return s + np.sin(2 * np.pi * np.cumsum(180 + 2600 * (t / d) ** 2) / SR) * (t / d) ** 1.8 * 0.22

    chords = [(55.0, [220.0, 261.63, 329.63]), (43.65, [174.61, 220.0, 261.63]),
              (65.41, [261.63, 329.63, 392.0]), (49.0, [196.0, 246.94, 293.66])]
    bars = int(round(TOTAL / BAR))
    # chaos intro (the first two bars): stuttering kicks that speed up, noise bursts, a long riser into the drop
    for t0 in (0.0, 0.5, 0.95, 1.3, 1.6, 1.85, 2.05, 2.22, 2.37, 2.5, 2.61, 2.71, 2.8, 2.88, 2.95, 3.02, 3.08, 3.13, 3.18, 3.22, 3.26):
        put(t0, kick(0.25), 0.9)
    for t0 in (0.12, 0.7, 1.05, 1.45, 1.9, 2.3, 2.6):
        put(t0, hat(open_=True) * 1.4, 0.8, 0.5)
    put(0.6, riser(SCENE_LEN - 0.6), 1.0)
    burst = filt(rng.normal(0, 1, int(0.06 * SR)), "high", 2500) * np.exp(-tt(0.06) * 50) * 0.5
    for k in range(10):
        put(0.2 + k * 0.33, burst, 0.6, 0.9)
    # from bar 3: the drop
    duck = np.ones(n)
    for bar in range(2, bars):
        t_bar = bar * BAR
        bass_f, triad = chords[(bar - 2) % 4]
        for b in range(4):
            tb = t_bar + b * BEAT
            put(tb, kick(), 1.0)
            put(tb + BEAT / 2, hat(), 0.9, 0.7)
            if b in (1, 3):
                put(tb, clap(), 0.9, 0.9)
            if b == 3:
                put(tb + BEAT * 0.75, hat(open_=True), 0.5, 0.8)
            i0 = int(tb * SR)
            m = min(int(BEAT * SR), n - i0)
            if m > 0:
                duck[i0 : i0 + m] = np.minimum(duck[i0 : i0 + m], 1 - 0.8 * np.exp(-np.arange(m) / SR / 0.11))
        for e in range(8):
            f = bass_f * (2 if e % 4 == 3 else 1)
            dd = BEAT / 2 * 0.9
            s = filt(saw(f, dd) * 0.7 + np.sin(2 * np.pi * f * tt(dd)) * 0.6, "low", 420) * np.exp(-tt(dd) * 5)
            put(t_bar + e * BEAT / 2, s * 0.5, 0.9, 0.9)
        if bar >= 3:
            for s16 in range(16):
                f = triad[(s16 * 2 + (s16 // 4)) % 3] * (2 if s16 % 4 == 2 else 1)
                tp = t_bar + s16 * BEAT / 4
                p = pluck(f)
                put(tp, p, 0.7 if s16 % 2 else 0.9, 0.9 if s16 % 2 else 0.7)
                put(tp + BEAT * 0.75, p, 0.18, 0.28)
    # pad under the drop, pumped by the kick
    pad_l = np.zeros(n)
    pad_r = np.zeros(n)
    for bar in range(2, bars):
        _, triad = chords[(bar - 2) % 4]
        d = BAR + 0.4
        sl = sum(saw(f, d, de) for f in triad for de in (-0.004, 0.0)) / 6
        sr_ = sum(saw(f, d, de) for f in triad for de in (0.003, 0.007)) / 6
        i0 = int(bar * BAR * SR)
        env = np.minimum(1, tt(d) * 6) * np.minimum(1, (d - tt(d)) * 3)
        m = min(len(sl), n - i0)
        pad_l[i0 : i0 + m] += (filt(sl, "low", 1500) * env)[:m]
        pad_r[i0 : i0 + m] += (filt(sr_, "low", 1500) * env)[:m]
    left += pad_l * duck * 0.5
    right += pad_r * duck * 0.5
    # a boom on every cut, with a riser and a clap roll leading into it
    for s in range(1, len(SCENES)):
        t0 = s * SCENE_LEN
        put(t0, boom(), 1.0)
        put(t0 - 0.9, riser(0.9), 0.7)
        for k in range(8):
            put(t0 - 0.5 + k * 0.0625, clap(0.12), 0.35 + k * 0.05)
    put(SCENE_LEN, boom(1.6), 0.5)
    out = np.tanh(np.stack([left, right], axis=1) * 0.9)
    out /= max(1e-6, np.abs(out).max()) / 0.89
    return (out * 32767).astype(np.int16)


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
           "-f", "rawvideo", "-pix_fmt", "bgra", "-s", f"{W}x{H}", "-r", str(FPS), "-i", "-",
           "-i", str(wav),
           "-c:v", "libx264", "-preset", "medium", "-crf", "20", "-pix_fmt", "yuv420p", "-profile:v", "high",
           "-g", str(FPS * 2), "-bf", "2", "-movflags", "+faststart",
           "-c:a", "aac", "-b:a", "192k", "-ar", str(SR), "-shortest", out]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    surf = cairo.ImageSurface(cairo.FORMAT_ARGB32, W, H)
    total = int(round(TOTAL * FPS))
    for f in range(total):
        t = f / FPS
        render(t, surf)
        buf = np.ndarray((H, W, 4), np.uint8, surf.get_data())
        lv = glitch_level(t)
        if lv > 0:
            rgb = glitch(np.ascontiguousarray(buf[:, :, [2, 1, 0]]), lv, f)
            frame = np.dstack([rgb[:, :, 2], rgb[:, :, 1], rgb[:, :, 0], np.full((H, W), 255, np.uint8)])
            proc.stdin.write(np.ascontiguousarray(frame).tobytes())
        else:
            proc.stdin.write(buf.tobytes())
    proc.stdin.close()
    proc.wait()
    wav.unlink()
    print(f"wrote {out}: {total} frames, {total / FPS:.1f}s, with sound")


def preview(times: list[float], path: str) -> None:
    from PIL import Image

    DATA.update(load_data())
    surf = cairo.ImageSurface(cairo.FORMAT_ARGB32, W, H)
    tiles = []
    for t in times:
        render(t, surf)
        buf = np.ndarray((H, W, 4), np.uint8, surf.get_data())
        rgb = np.ascontiguousarray(buf[:, :, [2, 1, 0]])
        if glitch_level(t) > 0:
            rgb = glitch(rgb, glitch_level(t), int(t * FPS))
        tiles.append(Image.fromarray(rgb).resize((640, 360)))
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
