# offset.offset_executive_p2: #855, the protective stitch closes the lobe the crossing-loop resolution dropped (2026-10-26)

Fork-only behaviour of the #854 opt-in (`Kernel::fixes.offset_resolve_crossing_loops`,
WR-345). Upstream at `c9c5979` does not resolve crossing loops, so there is no C++ result
to compare.

## Inputs

tangency-labs/solidcore#855 is a 40 x 40 x 10 plate with an r2 through hole at (x, 20),
shelled 1.5 inward with the top open. The hole's offset (r 3.5) crosses the side's offset
(x = 1.5) by `1.5 + 3.5 - x`. The probe is
`crates/solidcore/tests/issue_855_hole_offset_grazing_wall.rs`
(`issue_855_by_crossing_depth`, ignored).

## What happened

With the opt-in, the result depended on the crossing depth:

- 0.5 and 0.28: 15 faces, 36 edges, 24 vertices. A manifold solid, and the shell
  difference runs.
- 0.25, 0.1, 0.03, 0.01 and 0.001: 11 faces, 23 edges, 14 vertices.
  - Two vertical edges at the crossing lines have 8 edge uses each.
  - The offset is not a manifold solid, so `shell_brep` returns it raw (no difference).
  - The outer skin, the hole and the rim are missing.

## Cause

Traced with gdb (gdb command files `g855`, `g855b`, `g855c`).

- **Resolution and merge are the same at every depth.** `resolve_crossing_loops_of_face`
  resolves the floor and top offsets identically at 0.5 and 0.1: two crossings each, and
  the line stretch and the arc stretch between them dropped. The merged offset has the
  same 11 faces at both depths. Two of them are pieces the merge cut off:
  - the side offset's strip between the crossing lines (face 5);
  - the hole offset's sliver beyond x = 1.5 (face 6).
- **What differs is their regions.** At 0.5 both sides of both faces are the infinite
  region, and `MakeManifold` deletes them. At 0.1 a third region lies between them, so
  `MakeManifold` keeps both.
- **Where that region comes from.** The protective `StitchIntoSolid` after the merge
  (`offset_executive/p2.rs`, cpp:1772, a faithful port of `SmStitch.cpp:2659`) creates it:
  - `glue_edges_at_vertex` glues the strip's top and bottom edges to the sliver's;
  - the glue closes the lobe, and `test_closure_split_region` makes it a region of its own.

  The gap between those edges is the crossing depth. The stitch tolerance doubles from
  a tenth of the B-rep tolerance while it is below 1/200 of the bounding-box diagonal
  (about 53.1 here, so the last tolerance is about 0.262). Depths below that are glued and
  depths above are not. The depth sweep puts the switch between 0.28 and 0.25.

  Those edges are open on purpose: the floor and top no longer reach under the lobe,
  because the resolution dropped those stretches. They are not modelling gaps.

## Change (opt-in only)

1. `resolve_crossing_loops` records every stretch it drops: its end points and the
   middle of its parameter interval (`OffsetExecutive::crossing_open_stretches`).
2. Before the stitch, `edges_left_open` (`crossing_loops.rs`) finds the merged result's
   lamina edges along those stretches.
   - The merge rebuilds every edge, so ids cannot carry over, and no attribute
     distinguishes these pieces from the kept pieces beside them. The match is by points,
     using the model tolerances already set:
     - both end vertices within their own `Vertex::tolerance`, in either order;
     - the edge curve's middle within its own `Edge::tolerance`.
   - The match must be one to one. Otherwise no edge is kept open.
3. `stitch_into_solid_keeping_open` (`heal/stitch/p2.rs`, `pub(crate)`) runs the
   unchanged stitch, except that `glue_edges_at_vertex` never glues an edge of that set
   (`Stitch::keep_open`).
4. `MakeManifold` then removes the strip and the sliver, as it does at 0.5.

With the opt-in off, nothing is recorded and `stitch_into_solid` runs as before.

## Test

`issue_855_with_crossing_loops_resolved_by_depth` asserts, at crossing depths 0.5, 0.28,
0.25, 0.1 and 0.001, the 15-face shell: intact topology, a manifold solid, and the closed
form volume to 1e-6.

## Not done

The alternative was to trim the side and hole offsets against each other before the
merge, so the strip and sliver never exist. It was shelled: it duplicates the merge's
surface intersection.
