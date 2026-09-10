# Repaint a finished figure for a dark slide: knock the paper out to
# transparency, and recolour the ink so it reads light-on-dark.
#
# This is the stopgap path. The real one is fig_gydb_rt_layers.R, which now
# writes *_dark variants itself -- use that when ggtree is available. This
# script exists because the ink in a line figure is already separable from the
# paper by luminance alone, so a rendered PNG can be converted without R.
#
# How it works: the alpha channel is the ink, stretched so the darkest pixel in
# the figure becomes fully opaque. That is the point. A tree drawn in #8d8c87
# on white is only ~45% of the way to black, so a naive alpha = 1 - luminance
# leaves every branch half-transparent and the figure stays as faint on the
# slide as it was in the white box. Stretching to the observed ink floor makes
# the branches solid while keeping the anti-aliased edges soft.
#
#   python scripts/slide_transparent.py IN.png OUT.png [--ink "#f4f3f8"]
#                                       [--grow N] [--no-crop]
#
#   --grow N   dilate the ink by N pixels. Hairlines rendered at 400 dpi are
#              2-3 px wide and thin out when a slide scales them down; 1-2 is
#              usually enough. Check dense regions do not blob together.

import argparse

import numpy as np
from PIL import Image


def hex_to_rgb(s):
    s = s.lstrip("#")
    return tuple(int(s[i:i + 2], 16) for i in (0, 2, 4))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--ink", default="#f4f3f8")
    ap.add_argument("--grow", type=int, default=0)
    ap.add_argument("--no-crop", action="store_true")
    args = ap.parse_args()

    im = Image.open(args.src).convert("RGB")
    rgb = np.asarray(im).astype(np.float32)
    lum = rgb @ np.array([0.2126, 0.7152, 0.0722], dtype=np.float32)

    # Paper is whatever the bulk of the image is; ink floor is the darkest 0.1%
    # rather than the single darkest pixel, so one stray dark speck cannot set
    # the scale for the whole figure.
    paper = np.percentile(lum, 99.0)
    floor = np.percentile(lum, 0.1)
    if paper - floor < 1.0:
        raise SystemExit("no ink found: image is nearly uniform")

    alpha = np.clip((paper - lum) / (paper - floor), 0.0, 1.0)

    if args.grow > 0:
        from PIL import ImageFilter
        a = Image.fromarray((alpha * 255).astype(np.uint8))
        a = a.filter(ImageFilter.MaxFilter(2 * args.grow + 1))
        alpha = np.asarray(a).astype(np.float32) / 255.0

    h, w = alpha.shape
    out = np.zeros((h, w, 4), dtype=np.uint8)
    out[..., 0], out[..., 1], out[..., 2] = hex_to_rgb(args.ink)
    out[..., 3] = (alpha * 255).round().astype(np.uint8)

    img = Image.fromarray(out, mode="RGBA")

    if not args.no_crop:
        box = img.getchannel("A").point(lambda v: 255 if v > 2 else 0).getbbox()
        if box:
            img = img.crop(box)

    img.save(args.dst)
    print(f"wrote {args.dst}  {img.size[0]}x{img.size[1]}  "
          f"paper={paper:.0f} floor={floor:.0f} ink={args.ink}")


if __name__ == "__main__":
    main()
