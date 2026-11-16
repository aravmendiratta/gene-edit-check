# Design language

How gene-edit-check looks and sounds, so a figure, a terminal table and a README section feel like the same
product. The tokens live in code (`src/gene_edit_check/ui/plots.py`, `PALETTES`); this page explains the choices.

## Principles

1. **Quiet by default.** White or near-black surface, grey for everything that is context, one accent for the
   thing being discussed. If everything is coloured, nothing is.
2. **Colour has one job.** Blue is the subject. Orange means a check failed or a limit was crossed. Nothing else
   is orange. There is no green and no red.
3. **Shape before colour.** A pass is a tick and a fail is a cross, so the figures read in greyscale and for
   colour-blind readers.
4. **Say the finding.** Titles state the result ("No candidate passed every check"), subtitles say how to read
   the chart, and only the few points that matter get a label.
5. **Honest by construction.** Every figure footer says "model predictions, not lab measurements", and
   uncertainty is drawn when it exists.

## Palette

| Token | Light | Dark | Use |
|---|---|---|---|
| Surface | `#ffffff` | `#0d1117` | Page and figure background |
| Panel | `#f7f7f5` | `#151b23` | Cards and banded rows |
| Ink | `#111111` | `#f2f2f0` | Titles and values |
| Ink 2 | `#55534f` | `#b4b3ae` | Subtitles, axis labels |
| Ink 3 | `#8a8985` | `#8a8e96` | Footers, tick labels |
| Hairline | `#ebebe8` | `#262c35` | Grid and borders |
| Neutral | `#b9b8b3` | `#5d626b` | Context marks |
| Blue | `#2a78d6` | `#5aa0f2` | Accent: the subject |
| Orange | `#eb6834` | `#f08457` | Failure or limit |

The two accent hues are separated by a Delta E of 15 or more under colour-blind simulation (checked with the
dataviz palette validator). The grey is context, not a data series, so it is exempt from the lightness band.
The dark values are lighter steps of the same hues, not an inversion.

## Marks

- Thin lines (1 to 2 pt), hairline gridlines, no chart borders, no dual axes.
- Markers of at least 8 pt with a surface-coloured ring where they overlap.
- Log axes carry plain tick labels (0.01, 0.1, 1), never powers of ten.
- Ticks and crosses for pass and fail; filled circle for "all edits", open circle for "random edits only".

## Type

Inter, then SF Pro Display, Segoe UI, Helvetica Neue, Arial. Titles are semibold at 15 pt, subtitles regular at
about 10 pt in Ink 2, footers 7.6 pt in Ink 3. Spacing around titles is set in inches so it holds at any height.

## Brand

- **Mark:** a rounded blue square holding a white tick. `assets/logo.svg`.
- **Wordmark:** `gene-edit-` in Ink, `check` in Blue, semibold.
- **Banner and social card:** `assets/banner.svg` and `banner-dark.svg` and
  `assets/social-card.png` (1280 x 640, for the GitHub social preview).
- Regenerate with `python scripts/make_brand.py`.

## Light and dark

Every README figure is written twice, `fig_x.png` and `fig_x-dark.png`, and the README serves the right one
with a `<picture>` element. Switch palettes in code with `plots.use_mode("dark")`.

## Terminal

The CLI uses the same two hues: blue for the subject and passes, orange for failures and flags, dim grey for
everything else. Output is UTF-8 so the tables survive any console.

## Voice

Plain sentences, no jargon without a gloss, no exclamation marks, no claims the data do not support. A
negative result is stated as plainly as a positive one.
