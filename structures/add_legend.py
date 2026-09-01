"""
Post-process the AF3-vs-Boltz2 comparison renders to add a color legend
(blue = AF3, orange = Boltz-2), since PyMOL has no native legend widget.
Colors match exactly what was used in visualize_af3_vs_boltz_seeds.py
(PyMOL 'skyblue' and 'orange').

Usage: python3 structures/add_legend.py  (run with a python env that has Pillow)
"""
import os

from PIL import Image, ImageDraw, ImageFont

BASE = "/Users/ivanatang/Developer/SH2/structures"
IMAGES = [
    "af3_vs_boltz_seeds_full.png",
    "af3_vs_boltz_seeds_peptide_closeup.png",
    "af3_vs_boltz_seeds_ptyr_closeup.png",
]

AF3_COLOR = (51, 128, 204)     # pymol 'skyblue' (0.2, 0.5, 0.8) * 255
BOLTZ_COLOR = (255, 128, 0)    # pymol 'orange'  (1.0, 0.5, 0.0) * 255

for fname in IMAGES:
    path = f"{BASE}/{fname}"
    img = Image.open(path).convert("RGBA")
    draw = ImageDraw.Draw(img)

    # scale legend to image size
    w, h = img.size
    font_size = max(24, w // 45)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", font_size)
    except OSError:
        font = ImageFont.load_default()

    swatch = font_size
    pad = font_size // 2
    line_gap = int(font_size * 1.6)
    x0, y0 = pad, pad

    entries = [("AF3", AF3_COLOR), ("Boltz-2", BOLTZ_COLOR)]

    # measure legend box size
    text_widths = [draw.textlength(label, font=font) for label, _ in entries]
    box_w = int(swatch + pad * 2 + max(text_widths) + pad)
    box_h = int(pad + line_gap * len(entries))

    # semi-transparent white background box for legibility
    overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
    odraw = ImageDraw.Draw(overlay)
    odraw.rounded_rectangle(
        [x0 - pad // 2, y0 - pad // 2, x0 + box_w, y0 + box_h],
        radius=pad // 2,
        fill=(255, 255, 255, 220),
        outline=(120, 120, 120, 255),
        width=2,
    )
    img = Image.alpha_composite(img, overlay)
    draw = ImageDraw.Draw(img)

    for i, (label, color) in enumerate(entries):
        y = y0 + i * line_gap
        draw.rectangle([x0, y, x0 + swatch, y + swatch], fill=color + (255,), outline=(0, 0, 0, 255))
        draw.text((x0 + swatch + pad, y - 2), label, fill=(0, 0, 0, 255), font=font)

    img.convert("RGB").save(path)
    print(f"legend added -> {path}")
