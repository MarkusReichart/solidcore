"""VectorCraft steps for #2978: a plane beside a bare net, in one normal section of the edge.
The net's trace runs up from the edge point p to the side of the net's box, and on past it as
the extension the executive builds (Ruled: the line along the side's tangent). The moved
plane's trace x0 + tau u meets that line; Krawczyk's box around the meeting proves it, and the
reach is sqrt(|x0 - p|^2 + tau^2).
Usage: python make_vc_steps.py OUT.json
"""
import json
import math
import sys

BLUE, BLACK, PALE, GREY, DARK, FACE = "#1f6fb2", "#000000", "#dce9f5", "#999999", "#555555", "#f0f0f0"
FONT = "Helvetica"
W, H = 980, 560

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


def line(a, b):
    cmd("shape.line", x1=a[0], y1=a[1], x2=b[0], y2=b[1])


def polyline(points):
    cmd("path.create", anchors=[{"x": x, "y": y} for x, y in points], closed=False)


def dot(c, r=4.5, color=BLACK):
    style(fill=color)
    cmd("shape.ellipse", x=c[0] - r, y=c[1] - r, width=2 * r, height=2 * r)


def arrowhead(tip, frm, color=BLACK):
    dx, dy = tip[0] - frm[0], tip[1] - frm[1]
    n = math.hypot(dx, dy)
    ux, uy = dx / n, dy / n
    a = (tip[0] - 10 * ux + 5 * uy, tip[1] - 10 * uy - 5 * ux)
    b = (tip[0] - 10 * ux - 5 * uy, tip[1] - 10 * uy + 5 * ux)
    style(fill=color)
    cmd("path.create", anchors=[{"x": tip[0], "y": tip[1]}, {"x": a[0], "y": a[1]}, {"x": b[0], "y": b[1]}], closed=True)


cmd("file.new", width=W, height=H, units="Points")
# white background
style(fill="#ffffff")
cmd("shape.rectangle", x=0, y=0, width=W, height=H)
text(24, 44, "A plane beside a net, in one normal section of the edge", size=20, bold=True)

# the net's trace: from the edge point p up to the side of its box, a curve
p = (330.0, 470.0)
side = (420.0, 250.0)


def trace(t):
    # a quadratic from p to the side, leaning right
    c = (300.0, 340.0)
    return ((1 - t) ** 2 * p[0] + 2 * t * (1 - t) * c[0] + t * t * side[0],
            (1 - t) ** 2 * p[1] + 2 * t * (1 - t) * c[1] + t * t * side[1])


net = [trace(i / 40) for i in range(41)]
# the tangent at the side: 2 (side - c)
tx, ty = 2 * (side[0] - 300.0), 2 * (side[1] - 340.0)
tn = math.hypot(tx, ty)
tx, ty = tx / tn, ty / tn
ext_end = (side[0] + 215 * tx, side[1] + 215 * ty)

# the faces, pale, on the solid's side of the net's trace
style(fill=FACE)
cmd("path.create", anchors=[{"x": x, "y": y} for x, y in net] + [{"x": 170, "y": 250}, {"x": 170, "y": 470}], closed=True)
text(190, 440, "the solid", size=12, color=DARK)

style(stroke=BLUE, weight=3)
polyline(net)
style(stroke=BLUE, weight=3, dash=[9, 6])
line(side, ext_end)

# the side of the net's box
style(stroke=GREY, weight=1.5, dash=[6, 4])
line((side[0] - 34, side[1]), (side[0] + 34, side[1]))
dot(side, r=3.5, color=DARK)
text(side[0] + 14, side[1] + 22, "the side of the net's box", size=12)

# the moved plane's trace, a line
a = (150.0, 150.0)
b = (860.0, 118.0)
ux, uy = b[0] - a[0], b[1] - a[1]
un = math.hypot(ux, uy)
ux, uy = ux / un, uy / un
style(stroke=BLACK, weight=1.5)
line(a, b)
text(690, 100, "the moved plane's trace  x0 + τ u", size=14)

# the foot x0 of p on the line, and the meeting m with the extension
s = (p[0] - a[0]) * ux + (p[1] - a[1]) * uy
x0 = (a[0] + s * ux, a[1] + s * uy)
# solve a + t u = side + r (tx, ty)
det = ux * (-ty) - uy * (-tx)
t = ((side[0] - a[0]) * (-ty) - (side[1] - a[1]) * (-tx)) / det
m = (a[0] + t * ux, a[1] + t * uy)

style(stroke=GREY, weight=1.5, dash=[6, 4])
line(p, x0)
text(250, 300, "|x0 − p|", size=14)

# tau along the line from x0 to m, offset above it
off = (uy * 22, -ux * 22)
style(stroke=DARK, weight=1.25)
line((x0[0] + off[0], x0[1] + off[1]), (m[0] + off[0] - 4 * ux, m[1] + off[1] - 4 * uy))
arrowhead((m[0] + off[0], m[1] + off[1]), (x0[0] + off[0], x0[1] + off[1]), color=DARK)
text((x0[0] + m[0]) / 2 - 4, (x0[1] + m[1]) / 2 - 30, "τ", size=14)

# the Krawczyk box around the meeting
style(stroke=DARK, weight=1.0, fill=PALE)
cmd("shape.rectangle", x=m[0] - 18, y=m[1] - 13, width=36, height=26)
style(stroke=BLACK, weight=1.5)
line((m[0] - 18, m[1] - 18 * uy / ux), (m[0] + 18, m[1] + 18 * uy / ux))
style(stroke=BLUE, weight=3, dash=[9, 6])
line((m[0] - 13 * tx, m[1] - 13 * ty), (m[0] + 13 * tx, m[1] + 13 * ty))
dot(m, r=3.5)
style(stroke=GREY, weight=1.0)
line((m[0] + 20, m[1] + 12), (600, 210))
text(606, 206, "Krawczyk's box Y: K(Y) inside Y", size=14)
text(606, 226, "proves one meeting in it, for every", size=12)
text(606, 242, "section of a piece of the edge", size=12)

dot(p)
text(p[0] + 12, p[1] + 6, "p, a point of the edge", size=14)
dot(x0, r=3.5, color=DARK)
text(x0[0] - 30, x0[1] - 14, "x0", size=14)

text(352, 400, "the net's trace", size=14, color=BLACK)
text(330, 208, "past the side: the", size=12)
text(330, 224, "Ruled extension's line", size=12)

text(560, 420, "reach ≤ √(|x0 − p|² + τ²), the box's τ enclosed;", size=14)
text(560, 442, "|τ| is the arc along the moved plane's trace", size=12)
text(560, 462, "(over each piece of the edge's parameters,", size=12)
text(560, 478, "the foot and the direction by their mean values)", size=12)

json.dump(steps, open(sys.argv[1], "w", encoding="utf-8", newline="\n"), indent=0)
print(len(steps), "steps")
