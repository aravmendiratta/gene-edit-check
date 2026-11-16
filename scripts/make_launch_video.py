"""Build a 22-second launch video from the project's own figures. No audio, no network.

    pip install pillow imageio imageio-ffmpeg
    python scripts/make_launch_video.py            # writes assets/launch.mp4

Every claim on screen is one the README makes. The story is the honest one: no edit passed, and the
model misses the published enhancers.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from gene_edit_check.ui import brand

W, H, FPS = 1280, 720, 24
PAPER, INK, MUTE, TEAL, CORAL = (244, 238, 224), (27, 38, 34), (91, 101, 95), (31, 122, 104), (224, 80, 44)
OUT = Path("outputs")

_SERIF = ["georgia.ttf", "Georgia.ttf", "IowanOldStyle.ttf", "DejaVuSerif.ttf"]
_SANS = ["segoeui.ttf", "Inter-Regular.ttf", "Helvetica.ttc", "Arial.ttf", "DejaVuSans.ttf"]


def font(names: list[str], size: int) -> ImageFont.FreeTypeFont:
    for n in names:
        try:
            return ImageFont.truetype(n, size)
        except OSError:
            continue
    return ImageFont.load_default()


def icon_png(size: int) -> Image.Image:
    """Rasterise the helix mark through matplotlib so the video needs no SVG renderer."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(size / 100, size / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(100, 0)
    ax.axis("off")
    fig.patch.set_alpha(0)
    ax.patch.set_alpha(0)
    brand.draw_mark(ax, 0, 0, 100)
    fig.canvas.draw()
    img = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).convert("RGBA")
    plt.close(fig)
    # the corners outside the rounded tile are transparent in the figure; paste on paper later
    return img


def fit(img: Image.Image, max_w: int, max_h: int) -> Image.Image:
    r = min(max_w / img.width, max_h / img.height)
    return img.resize((int(img.width * r), int(img.height * r)), Image.LANCZOS)


def ease(t: float) -> float:
    t = max(0.0, min(1.0, t))
    return t * t * (3 - 2 * t)


def text_center(d: ImageDraw.ImageDraw, y: int, s: str, f: ImageFont.FreeTypeFont, fill: tuple, alpha: float = 1.0) -> None:
    w = d.textlength(s, font=f)
    d.text(((W - w) / 2, y), s, font=f, fill=fill + (int(255 * alpha),))


def scene(i: int, t: float, dur: float, assets: dict) -> Image.Image:
    """Render scene `i` at local time `t` seconds onto a paper background."""
    base = Image.new("RGBA", (W, H), PAPER + (255,))
    layer = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    d = ImageDraw.Draw(layer)
    a_in = ease(t / 0.6)
    a_out = 1 - ease((t - (dur - 0.5)) / 0.5)
    a = min(a_in, a_out)
    serif, sans = assets["serif"], assets["sans"]
    if i == 0:
        s = 0.86 + 0.14 * ease(t / 1.0)
        ic = assets["icon"].resize((int(220 * s), int(220 * s)), Image.LANCZOS)
        layer.alpha_composite(ic, ((W - ic.width) // 2, 150 + int((220 - ic.height) / 2)))
        f = assets["serif_xl"]
        w1 = d.textlength("gene-edit-", font=f)
        w2 = d.textlength("check", font=f)
        x0 = (W - (w1 + w2)) / 2
        d.text((x0, 420), "gene-edit-", font=f, fill=INK + (255,))
        d.text((x0 + w1, 420), "check", font=f, fill=TEAL + (255,))
        text_center(d, 548, "Check a gene edit on a computer before the lab.", assets["sans_m"], MUTE, ease((t - 0.8) / 0.6))
    elif i == 1:
        text_center(d, 230, "Testing a gene edit in a lab takes weeks.", assets["serif_l"], INK, ease(t / 0.6))
        text_center(d, 330, "So check it on a computer first.", assets["serif_l"], TEAL, ease((t - 1.6) / 0.6))
    elif i == 2:
        text_center(d, 70, "Propose on some patients. Judge on others.", assets["serif_m"], INK)
        im = fit(assets["pipeline"], 1160, 420)
        layer.alpha_composite(im, ((W - im.width) // 2, 220))
    elif i == 3:
        text_center(d, 40, "8 edits. 6 safety checks. None passed.", assets["serif_m"], INK)
        text_center(d, 100, "7 of 8 disturbed DNA methylation.", sans, CORAL, ease((t - 1.2) / 0.5))
        im = fit(assets["scorecard"], 1000, 540)
        layer.alpha_composite(im, ((W - im.width) // 2, 160))
    elif i == 4:
        text_center(d, 60, "We checked the model against a lab screen.", assets["serif_m"], INK)
        text_center(d, 120, "It finds elements near the promoter, not the +15 kb enhancers.", sans, MUTE, ease((t - 0.8) / 0.5))
        im = fit(assets["enhancers"], 1160, 480)
        layer.alpha_composite(im, ((W - im.width) // 2, 195))
    else:
        ic = assets["icon"].resize((120, 120), Image.LANCZOS)
        layer.alpha_composite(ic, ((W - 120) // 2, 130))
        text_center(d, 290, "pip install gene-edit-check", assets["mono"], INK)
        text_center(d, 380, "github.com/aravmendiratta/gene-edit-check", sans, TEAL)
        text_center(d, 560, "Research use only. Model predictions, not lab measurements.", assets["sans_s"], MUTE)
    if a < 1:
        r, g, b, al = layer.split()
        layer.putalpha(al.point(lambda v: int(v * a)))
    base.alpha_composite(layer)
    return base.convert("RGB")


def main(out: str = "assets/launch.mp4") -> None:
    import imageio.v2 as imageio

    assets = {
        "serif": font(_SERIF, 40), "serif_m": font(_SERIF, 44), "serif_l": font(_SERIF, 60), "serif_xl": font(_SERIF, 92),
        "sans": font(_SANS, 30), "sans_m": font(_SANS, 32), "sans_s": font(_SANS, 22),
        "mono": font(["consola.ttf", "Menlo.ttc", "DejaVuSansMono.ttf"], 48),
        "icon": icon_png(440),
        "pipeline": Image.open(OUT / "fig_pipeline.png").convert("RGBA"),
        "scorecard": Image.open(OUT / "fig_scorecard.png").convert("RGBA"),
        "enhancers": Image.open(OUT / "fig_enhancers.png").convert("RGBA"),
    }
    # the figures sit on a slightly lighter card so they read against the paper background
    plan = [(0, 3.4), (1, 4.2), (2, 3.8), (3, 5.2), (4, 5.0), (5, 3.4)]
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    writer = imageio.get_writer(out, fps=FPS, codec="libx264", quality=8, macro_block_size=1,
                                ffmpeg_params=["-movflags", "+faststart"])
    n = 0
    for i, dur in plan:
        for k in range(int(dur * FPS)):
            writer.append_data(np.asarray(scene(i, k / FPS, dur, assets)))
            n += 1
    writer.close()
    print(f"wrote {out}: {n} frames, {n / FPS:.1f}s")


if __name__ == "__main__":
    main(*sys.argv[1:2])
