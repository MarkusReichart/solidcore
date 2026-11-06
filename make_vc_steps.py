"""VectorCraft steps for #2872: where the split at a radius law's corners runs.
Before: only FilletSession::do_fillet split the round; solidcore-ops drives the session step
by step (create_corners, its corner settings, do_filleting) and its round stayed one face that
kinks at the corner. After: do_filleting ends with the split, so both drivers reach it.
Usage: python make_vc_steps.py OUT.json
"""
import json
import sys

BLUE, BLACK, PALE, GREY, DARK, FACE = "#1f6fb2", "#000000", "#dce9f5", "#999999", "#555555", "#f0f0f0"
FONT = "Helvetica"
W, H = 980, 560
BW, BH = 128, 46

steps = []


def cmd(c, **p):
    steps.append({"command": c, "params": p})


def style(stroke=None, weight=1.0, dash=None, fill=None):
    cmd("select.none")
    if fill:
        cmd("paint.setFill", color=fill)
    else:
        cmd("paint.setFill", none=True)
    if stroke:
        cmd("paint.setStroke", color=stroke)
        cmd("stroke.set", weight=weight, dash=dash, cap="round", join="round")
    else:
        cmd("paint.setStroke", none=True)


def text(x, y, s, size=14, bold=False, color=BLACK):
    cmd("select.none")
    p = {"x": x, "y": y, "text": s, "size": size, "font": FONT, "color": color, "justify": "left"}
    if bold:
        p["style"] = "Bold"
    cmd("text.create", **p)


def box(x, y, lines, subject=False, missing=False):
    if missing:
        style(stroke=GREY, weight=1.5, dash=[6, 4], fill="#ffffff")
    elif subject:
        style(stroke=BLUE, weight=3, fill=PALE)
    else:
        style(stroke=DARK, weight=1.0, fill=FACE)
    cmd("shape.rectangle", x=x, y=y, width=BW, height=BH)
    top = y + BH / 2 - 7 * (len(lines) - 1) + 4
    for k, line in enumerate(lines):
        text(x + 8, top + 14 * k, line, size=12, color=GREY if missing else BLACK)


def arrow(x1, y1, x2, y2):
    style(stroke=BLACK, weight=1.25)
    cmd("shape.line", x1=x1, y1=y1, x2=x2 - 7, y2=y2)
    style(fill=BLACK)
    cmd("path.create", anchors=[{"x": x2, "y": y2}, {"x": x2 - 9, "y": y2 - 5}, {"x": x2 - 9, "y": y2 + 5}], closed=True)


def lane(x0, y, title, boxes):
    text(x0, y - 12, title, size=14)
    x = x0
    for k, b in enumerate(boxes):
        box(x, y, **b)
        if k + 1 < len(boxes):
            arrow(x + BW, y + BH / 2, x + BW + 30, y + BH / 2)
        x += BW + 30


cmd("file.new", width=W, height=H, units="Points")
style(fill="#ffffff")
cmd("shape.rectangle", x=0, y=0, width=W, height=H)

# before
L = 24
text(L, 44, "Before", size=20, bold=True)
lane(L, 110, "FilletSession::do_fillet", [
    {"lines": ["create_corners"]},
    {"lines": ["do_filleting"]},
    {"lines": ["split at the", "law's corners"], "subject": True},
])
lane(L, 250, "solidcore-ops fillet_variable (step by step)", [
    {"lines": ["create_corners"]},
    {"lines": ["corner settings"]},
    {"lines": ["do_filleting"]},
])
box(L + 2 * (BW + 30), 330, ["no split"], missing=True)
style(stroke=GREY, weight=1.5, dash=[6, 4])
cmd("shape.line", x1=L + 2 * (BW + 30) + BW / 2, y1=296, x2=L + 2 * (BW + 30) + BW / 2, y2=330)
text(L, 420, "The round stays one face that kinks at the law's corner:", size=12)
text(L, 438, "not a valid solid, refused NotSolid. Before that, ops", size=12)
text(L, 456, "refused the law outright (Unsupported).", size=12)

# the divider
style(stroke=GREY, weight=1.0)
cmd("shape.line", x1=W / 2, y1=30, x2=W / 2, y2=H - 30)

# after
R = W / 2 + 24
text(R, 44, "After", size=20, bold=True)
lane(R, 110, "FilletSession::do_fillet", [
    {"lines": ["create_corners"]},
    {"lines": ["do_filleting"]},
])
lane(R, 250, "solidcore-ops fillet_variable (step by step)", [
    {"lines": ["create_corners"]},
    {"lines": ["corner settings"]},
    {"lines": ["do_filleting"]},
])
# do_filleting expanded, below
ex = R
ey = 370
text(ex, ey - 12, "do_filleting", size=14)
box(ex, ey, ["filleting"])
arrow(ex + BW, ey + BH / 2, ex + BW + 30, ey + BH / 2)
box(ex + BW + 30, ey, ["split at the", "law's corners"], subject=True)
arrow(ex + 2 * BW + 30, ey + BH / 2, ex + 2 * BW + 60, ey + BH / 2)
box(ex + 2 * (BW + 30), ey, ["blend face log"])
# both drivers' do_filleting lead to it
style(stroke=GREY, weight=1.5, dash=[6, 4])
cmd("shape.line", x1=R + (BW + 30) + BW / 2, y1=156, x2=R + (BW + 30) + BW / 2, y2=200)
cmd("shape.line", x1=R + 2 * (BW + 30) + BW / 2, y1=296, x2=R + 2 * (BW + 30) + BW / 2, y2=352)
text(R + (BW + 30) + BW / 2 + 6, 190, "the same do_filleting", size=12)
text(R, 456, "Both drivers split the round at the law's corners:", size=12)
text(R, 474, "the peak 2, 4, 2 is 8 faces, removing 39.540969.", size=12)

json.dump(steps, open(sys.argv[1], "w", encoding="utf-8", newline="\n"), indent=0)
print(len(steps), "steps")
