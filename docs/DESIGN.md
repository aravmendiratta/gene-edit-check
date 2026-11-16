# Design language

How gene-edit-check looks and sounds, so a figure, a terminal table and a README section feel like the same
product. The tokens live in code (`src/gene_edit_check/ui/plots.py`, `PALETTES`); this page explains the choices.

## Principles

1. **Feel drawn, not generated.** A warm paper surface, a serif for titles and a mark with a pulse of life in
   it. The product is biology, so it should look closer to a field guide than to a dashboard.
2. **One accent, one warning.** Deep teal is the subject. Coral means a check failed or a limit was crossed,
   and it is also the colour of the edited base in the logo. Nothing else is coral. There is no green and no red.
3. **Shape before colour.** A pass is a tick and a fail is a cross, so figures read in greyscale and for
   colour-blind readers.
4. **Say the finding.** Titles state the result ("No candidate passed every check"), subtitles say how to read
   the chart, and only the few points that matter get a label.
5. **Honest by construction.** Every figure footer says "model predictions, not lab measurements", and
   uncertainty is drawn when it exists.

## The mark

A double helix, leaning and half-turned, with one rung in coral: the base being edited. The ribbons swell and
thin as each strand turns toward the viewer, and the colour drifts from cream to mint along the length, so it
reads as a single drawn object and not a flat icon. It sits on a deep forest tile for use as an app icon and
as a small stamp in the corner of every figure.

- Built in code, once, in `src/gene_edit_check/ui/brand.py`; drawn as SVG (`assets/logo.svg`, the banners) and
  as matplotlib patches (figures, social card).
- The wordmark is set in a serif: `gene-edit-` upright in ink, `check` in italic teal.

## Palette

| Token | Light | Dark | Use |
|---|---|---|---|
| Surface | `#fbf8f1` | `#0e1714` | Page and figure background (warm paper, deep forest) |
| Panel | `#f3eee2` | `#15211d` | Cards and banded rows |
| Ink | `#1b2622` | `#efe9dc` | Titles and values |
| Ink 2 | `#505a54` | `#b5b8ad` | Subtitles, axis labels |
| Ink 3 | `#84897f` | `#868c83` | Footers, tick labels |
| Hairline | `#e8e2d3` | `#24322d` | Grid and borders |
| Neutral | `#bdb6a6` | `#5b665f` | Context marks |
| Teal | `#1f7a68` | `#4fbf9f` | Accent: the subject |
| Coral | `#e0502c` | `#ff7a4d` | Failure or limit |
| Forest | `#12453b` | same | Icon tile |

The two accent hues are teal and coral; the grey is context, not a data series. The dark values are lighter
steps of the same hues, not an inversion.

## Marks

- Thin lines (1 to 2 pt), hairline gridlines, no chart borders, no dual axes.
- Markers of at least 8 pt with a surface-coloured ring where they overlap.
- Log axes carry plain tick labels (0.01, 0.1, 1), never powers of ten.
- Ticks and crosses for pass and fail; filled circle for "all edits", open circle for "random edits only".

## Type

Titles: Iowan Old Style, Source Serif 4, Palatino, Georgia, then the system serif, regular at 16.5 pt.
Everything else: Inter, SF Pro, Segoe UI, Helvetica Neue, Arial. Subtitles about 10 pt in Ink 2, footers
7.6 pt in Ink 3. Spacing around titles is set in inches so it holds at any figure height.

## Brand files

`assets/logo.svg` (icon), `assets/banner.svg` and `banner-dark.svg`, `assets/social-card.png` (1280 x 640, upload
it under repository Settings, Social preview). Regenerate with `python scripts/make_brand.py`.

## Light and dark

Every README figure is written twice, `fig_x.png` and `fig_x-dark.png`, and the README serves the right one
with a `<picture>` element. Switch palettes in code with `plots.use_mode("dark")`.

## Terminal

The CLI uses the same two hues: blue for the subject and passes, orange for failures and flags, dim grey for
everything else. Output is UTF-8 so the tables survive any console.

## Voice

Plain sentences, no jargon without a gloss, no exclamation marks, no claims the data do not support. A
negative result is stated as plainly as a positive one.
