#!/usr/bin/env python3
# Render boolean2 arrangement-trace phases to SVG (rasterized to PNG by chrome).
# Pure stdlib. viewBox is explicit per call, so "zoom all the way in" is just a
# small half-width.
#   python3 render_trace.py TRACE.json OUTDIR CX CY HALF1 [HALF2 ...]
# Optional OUTDIR/highlights.json: [{"x":..,"y":..,"label":..,"color":..}] rings
# specific verts (e.g. the imbalanced pair) on every phase.
import json
import os
import sys

W = 1100
MARGIN = 0.06

# color, stroke-width, directed(arrowhead), dashed  -- keyed on segment "kind"
KIND = {
    "input_edge": ("#c2c2c2", 1.6, False, False),
    "collapsed_edge": ("#9a9a9a", 1.6, False, False),
    "candidate_edge": ("#e6c6c6", 1.0, False, True),
    "listed_subsegment_unsplit": ("#b4b4b4", 1.4, False, False),
    "arrangement_subsegment": ("#4a90d9", 2.2, False, False),
    "arrangement_subsegment_unsplit": ("#aaccef", 1.2, False, True),
    "canonical_subedge": ("#1f9e89", 2.2, False, False),
    "retained_edge": ("#1a9850", 3.0, True, False),
}


def seg_style(kind):
    return KIND.get(kind, ("#bbbbbb", 1.5, False, False))


def pt_style(kind):
    return {
        "output_vertex": ("#157a3c", 4.2),
        "arrangement_vertex": ("#333333", 3.8),
        "canonical_vertex": ("#1f7a6a", 3.6),
        "merged_vertex": ("#333333", 3.4),
        "list_vertex": ("#444444", 3.2),
        "input_vertex": ("#222222", 3.4),
    }.get(kind, ("#333333", 3.2))


def collect_bbox(phase):
    xs, ys = [], []
    for p in phase.get("points", []):
        xs.append(p["xy"][0]); ys.append(p["xy"][1])
    for s in phase.get("segments", []):
        xs += [s["a"][0], s["b"][0]]; ys += [s["a"][1], s["b"][1]]
    for poly in phase.get("polygons", []):
        for v in poly["verts"]:
            xs.append(v[0]); ys.append(v[1])
    return xs, ys


def render(phase, eps, cx, cy, half, title, highlights, label_hl):
    x0, x1, y0, y1 = cx - half, cx + half, cy - half, cy + half

    def X(x):
        return (x - x0) / (x1 - x0) * W

    def Y(y):
        return W - (y - y0) / (y1 - y0) * W

    o = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{W}">',
         f'<rect width="{W}" height="{W}" fill="white"/>',
         '<defs><marker id="ah" markerWidth="10" markerHeight="10" refX="8" '
         'refY="3" orient="auto"><path d="M0,0 L8,3 L0,6 Z" fill="#1a9850"/>'
         '</marker></defs>']

    for poly in phase.get("polygons", []):
        pts = " ".join(f"{X(v[0]):.2f},{Y(v[1]):.2f}" for v in poly["verts"])
        o.append(f'<polygon points="{pts}" fill="#ff7f0e" fill-opacity="0.12" '
                 f'stroke="#ff7f0e" stroke-width="1.6"/>')

    for s in phase.get("segments", []):
        color, wdt, directed, dashed = seg_style(s.get("kind", ""))
        ax, ay, bx, by = X(s["a"][0]), Y(s["a"][1]), X(s["b"][0]), Y(s["b"][1])
        mk = ' marker-end="url(#ah)"' if directed else ""
        da = ' stroke-dasharray="4,3"' if dashed else ""
        o.append(f'<line x1="{ax:.2f}" y1="{ay:.2f}" x2="{bx:.2f}" y2="{by:.2f}" '
                 f'stroke="{color}" stroke-width="{wdt}"{mk}{da}/>')

    show_lbl = half < 1.0
    for p in phase.get("points", []):
        color, r = pt_style(p.get("kind", ""))
        px, py = X(p["xy"][0]), Y(p["xy"][1])
        if -30 < px < W + 30 and -30 < py < W + 30:
            o.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="{r}" fill="{color}" '
                     f'stroke="white" stroke-width="0.7"/>')
            lbl = p.get("id", "")
            if show_lbl and lbl:
                o.append(f'<text x="{px+6:.1f}" y="{py-6:.1f}" font-size="13" '
                         f'font-family="monospace" fill="#555">{lbl}</text>')

    for h in highlights:
        px, py = X(h["x"]), Y(h["y"])
        if -30 < px < W + 30 and -30 < py < W + 30:
            o.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="13" fill="none" '
                     f'stroke="{h.get("color","#d00")}" stroke-width="2.5"/>')
            if label_hl:
                o.append(f'<text x="{px+16:.1f}" y="{py+5:.1f}" font-size="16" '
                         f'font-family="monospace" fill="{h.get("color","#d00")}">'
                         f'{h.get("label","")}</text>')

    eps_px = eps / (x1 - x0) * W
    bar_px, bar_lbl = ((120, f"eps={eps:.3g} (={eps_px:.2g}px, sub-pixel)")
                       if eps_px < 2 else (eps_px, f"eps={eps:.3g}"))
    yb = W - 38
    o.append(f'<line x1="40" y1="{yb}" x2="{40+bar_px:.1f}" y2="{yb}" '
             f'stroke="#000" stroke-width="3"/>')
    o.append(f'<text x="40" y="{yb-9}" font-size="15" font-family="monospace" '
             f'fill="#000">{bar_lbl}</text>')
    o.append(f'<text x="20" y="30" font-size="20" font-family="monospace" '
             f'fill="#000">{title}</text>')
    o.append(f'<text x="20" y="54" font-size="14" font-family="monospace" '
             f'fill="#666">window={x1-x0:.3g} @ ({cx:.6g},{cy:.6g})</text>')
    o.append("</svg>")
    return "\n".join(o)


def main():
    trace = json.load(open(sys.argv[1]))
    outdir = sys.argv[2]
    cx, cy = float(sys.argv[3]), float(sys.argv[4])
    halves = [float(h) for h in sys.argv[5:]]
    eps = trace.get("eps", 0.0)
    hl_path = os.path.join(outdir, "highlights.json")
    highlights = json.load(open(hl_path)) if os.path.exists(hl_path) else []
    for i, phase in enumerate(trace["phases"]):
        name = phase["name"]
        xs, ys = collect_bbox(phase)
        lbl = (name == "filtered_output_edges")
        if xs:
            bcx, bcy = (min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2
            bhalf = max(max(xs) - min(xs), max(ys) - min(ys), 1e-9) / 2 * (1 + MARGIN * 4)
            open(f"{outdir}/{i:02d}_{name}_overview.svg", "w").write(
                render(phase, eps, bcx, bcy, bhalf, f"{i}:{name} [overview]", highlights, lbl))
        for h in halves:
            tag = f"z{h:g}".replace(".", "p").replace("-", "m")
            open(f"{outdir}/{i:02d}_{name}_{tag}.svg", "w").write(
                render(phase, eps, cx, cy, h, f"{i}:{name} [{tag}]", highlights, lbl))
    print("rendered", len(trace["phases"]), "phases x", 1 + len(halves), "views")


if __name__ == "__main__":
    main()
