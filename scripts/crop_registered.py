# Crop a set of transparent PNGs to ONE shared box.
#
# The overlay slides come off ggsave on a canvas sized for the widest of them,
# so each one carries a wide transparent margin and the tree lands at maybe
# half the frame. Cropping is free resolution on the slide -- but only if every
# figure in the set is cropped to the SAME box. Crop each to its own ink and
# the trees no longer line up, which is the one thing these files exist to do.
#
# So: union the ink boxes, pad, apply that single box to all of them. Files are
# rewritten in place and must already be RGBA with a real alpha channel.
#
#   python scripts/crop_registered.py OUT_DIR FILE.png [FILE.png ...] [--pad N]

import argparse
import os

from PIL import Image


def ink_box(img, thresh=2):
    return img.getchannel("A").point(lambda v: 255 if v > thresh else 0).getbbox()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--pad", type=int, default=24)
    args = ap.parse_args()

    imgs = []
    for f in args.files:
        im = Image.open(f)
        if im.mode != "RGBA":
            raise SystemExit(f"{f} is {im.mode}, not RGBA -- no alpha to crop on")
        imgs.append((f, im))

    sizes = {im.size for _, im in imgs}
    if len(sizes) != 1:
        raise SystemExit(f"inputs are not the same canvas: {sizes}")

    boxes = [ink_box(im) for _, im in imgs]
    if any(b is None for b in boxes):
        raise SystemExit("one input is fully transparent")

    w, h = imgs[0][1].size
    box = (max(0, min(b[0] for b in boxes) - args.pad),
           max(0, min(b[1] for b in boxes) - args.pad),
           min(w, max(b[2] for b in boxes) + args.pad),
           min(h, max(b[3] for b in boxes) + args.pad))

    for f, im in imgs:
        im.crop(box).save(f)
        print(f"{os.path.basename(f):36s} -> {box[2]-box[0]}x{box[3]-box[1]}")
    print(f"shared box {box} from canvas {w}x{h}")


if __name__ == "__main__":
    main()
