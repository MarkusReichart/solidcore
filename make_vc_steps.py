"""VectorCraft steps for the #2441 schematic: the ball R 10 cut at z = 4 seen from above.
Its rim is two edges joined at (+-sqrt 84, 0, 4). Left: the case's one point selects the
edge y > 0, so the round covers half the rim and ends at both vertices. Right: a point on
each edge, the whole rim is rounded. The pale band is the round in plan, from its contact
on the cut face (radius sqrt 60) to its contact on the ball (radius 10 sqrt 60 / 8).
Usage: python make_vc_steps.py OUT.json
"""
import json
import math
import sys

BLUE, BLACK, PALE, GREY, DARK, FACE = "#1f6fb2", "#000000", "#dce9f5", "#999999", "#555555", "#f0f0f0"
FONT = "Helvetica"
W, H = 960, 560
S = 15.0  # points per mm
RIM = math.sqrt(84.0) * S
INNER = math.sqrt(60.0) * S
OUTER = 10.0 * math.sqrt(60.0) / 8.0 * S
BALL = 10.0 * S

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


def at(cx, cy, deg, r):
    a = math.radians(deg)
    return cx + r * math.cos(a), cy - r * math.sin(a)


def arc(cx, cy, d0, d1, r, n=60):
    return [{"x": x, "y": y} for x, y in (at(cx, cy, d0 + (d1 - d0) * i / n, r) for i in range(n + 1))]


def text(x, y, s, size=14, bold=False):
    cmd("select.none")
    p = {"x": x, "y": y, "text": s, "size": size, "font": FONT, "color": BLACK, "justify": "left"}
    if bold:
        p["style"] = "Bold"
    cmd("text.create", **p)


def dot(x, y, r=4.0, color=BLACK):
    style(fill=color)
    cmd("shape.ellipse", x=x - r, y=y - r, width=2 * r, height=2 * r)


PANELS = [
    {"cx": 240, "title": "One point: half the rim", "picked": [90.0], "band": (0.0, 180.0),
     "note": ["Edge y > 0 only: the round ends at both vertices.",
              "5 faces, volume 3272.764885312"]},
    {"cx": 720, "title": "A point on each edge: the whole rim", "picked": [90.0, 270.0], "band": (0.0, 360.0),
     "note": ["Both edges: the round goes all round.",
              "4 faces, volume 3261.518763866 (closed form 3261.518762521)"]},
]
CY = 280

cmd("file.new", width=W, height=H, units="Points")
style(fill="#ffffff")
cmd("shape.rectangle", x=0, y=0, width=W, height=H)

for p in PANELS:
    cx = p["cx"]
    # the ball's outline in plan
    style(fill=FACE)
    cmd("shape.ellipse", x=cx - BALL, y=CY - BALL, width=2 * BALL, height=2 * BALL)
    # the round in plan: an annulus sector between its two contacts
    d0, d1 = p["band"]
    style(fill=PALE)
    if d1 - d0 >= 360.0:
        cmd("path.create", anchors=arc(cx, CY, 0, 180, OUTER) + arc(cx, CY, 180, 0, INNER), closed=True)
        style(fill=PALE)
        cmd("path.create", anchors=arc(cx, CY, 180, 360, OUTER) + arc(cx, CY, 360, 180, INNER), closed=True)
    else:
        cmd("path.create", anchors=arc(cx, CY, d0, d1, OUTER) + arc(cx, CY, d1, d0, INNER), closed=True)
    # the axes, construction
    style(stroke=GREY, weight=1.5, dash=[6, 4])
    cmd("shape.line", x1=cx - BALL - 20, y1=CY, x2=cx + BALL + 20, y2=CY)
    cmd("shape.line", x1=cx, y1=CY + BALL + 14, x2=cx, y2=CY - BALL - 14)
    text(cx + BALL + 22, CY + 5, "x", size=12)
    text(cx + 6, CY - BALL - 16, "y", size=12)
    # the rim's two edges: selected ones blue and thick, the other dark grey and dashed
    for e0, e1 in ((0.0, 180.0), (180.0, 360.0)):
        mid = (e0 + e1) / 2
        if any(e0 < q < e1 for q in p["picked"]):
            style(stroke=BLUE, weight=3)
        else:
            style(stroke=DARK, weight=2, dash=[5, 4])
        cmd("path.create", anchors=arc(cx, CY, e0, e1, RIM))
    # the two vertices joining them
    for d in (0.0, 180.0):
        x, y = at(cx, CY, d, RIM)
        dot(x, y, 5.0)
    # the picked points
    for q in p["picked"]:
        x, y = at(cx, CY, q, RIM)
        style(stroke=BLACK, weight=1.5, fill="#ffffff")
        cmd("shape.ellipse", x=x - 6, y=y - 6, width=12, height=12)
    text(cx - 215, 40, p["title"], size=20, bold=True)
    for k, line in enumerate(p["note"]):
        text(cx - 215, H - 52 + 18 * k, line, size=12)

# labels on the left panel only, beside what they name
cx = PANELS[0]["cx"]
x, y = at(cx, CY, 180.0, RIM)
text(x + 12, y + 24, "(-sqrt 84, 0, 4)", size=12)
x, y = at(cx, CY, 0.0, RIM)
text(x - 100, y + 24, "(sqrt 84, 0, 4)", size=12)
x, y = at(cx, CY, 90.0, RIM)
text(x + 12, y - 12, "picked point", size=12)
bx, by = at(cx, CY, 140.0, (INNER + OUTER) / 2)
style(stroke=BLACK, weight=0.75)
cmd("shape.line", x1=bx, y1=by, x2=bx - 40, y2=by - 50)
text(bx - 120, by - 56, "round, in plan", size=12)
x, y = at(cx, CY, 305.0, BALL)
text(x + 6, y + 18, "edge y < 0, not picked", size=12)
x, y = at(cx, CY, 60.0, RIM)
text(x + 12, y - 8, "edge y > 0", size=12)

json.dump(steps, open(sys.argv[1], "w", encoding="utf-8", newline="\n"), indent=0)
print(len(steps), "steps")
