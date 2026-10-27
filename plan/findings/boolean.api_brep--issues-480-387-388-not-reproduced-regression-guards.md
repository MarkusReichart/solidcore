# boolean.api_brep: issues #480, #387 and #388 are not reproduced by the port; tests guard them (2026-10-26)

Three SolidCraft boolean reports (from its `kernel-smlib` rig and kernel tests) that the Rust
port does not reproduce. None has been reproduced on this harness's C++ reference, and none
points to a defect location in upstream, so they are not upstream issues. Each has a test
that guards the port's current, correct behaviour.

## #480: union of a ring and a lip with cylinder walls 1e-4 to 1e-3 apart

- **Inputs:** the base is an annulus r 3..5 about (0, 10), z 0..2. The lip is a half annulus
  r 4..5 about (e, 10), y >= 10, extruded over z 2..20, with its cap on the base's top.
- **Reported (SolidCraft):** at e = 1e-4 and 1e-3 the union fails ("not a solid"). At e = 0
  and 0.01 it is valid with V = 355.000 (analytic 354.99997). A rig case with two lips
  returned 1066.0 mm³ instead of 769.47, with no error.
- **Rust:** a valid union with V = 355.000 (to 1e-4 relative) at e = 0, 1e-4, 1e-3 and 0.01.
- **Test:** `crates/solidcore/tests/issue_480_ring_lip_union.rs`,
  `issue_480_ring_lip_union_at_every_offset`.

## #387: a box minus an exact rational-NURBS cylinder

- **Inputs:** the box (-14, 0, -14)..(14, 20, 14) minus a cylinder (base (0, 50, -60), axis Z,
  r 40, h 120) given as exact rational NURBS (two rational quadratic half-cylinder faces,
  4 spans each).
- **Reported (SolidCraft):** an invalid result. Edge and vertex tolerances reach 0.78 mm
  against a body tolerance of 0.022 (tests #59 and #82), and `ValidatePointers` fails. The
  same cut with the analytic cylinder is valid.
- **Rust:** with the cylinder made NURBS by `SmBrep::TurnToNURBS`
  (`topology::brep::p3::turn_to_nurbs`), the cut is a manifold solid with the analytic volume
  (to 1e-5 relative), as is the analytic cut. The NURBS form differs from the issue's
  flat-array import, so this does not prove the issue's exact input is fine.
- **Test:** `crates/solidcore/tests/issue_387_nurbs_cylinder_cut.rs`,
  `issue_387_box_minus_nurbs_cylinder_matches_the_analytic_cut`.

## #388: an island inside a void, then cut

- **Inputs:** a 20 cube joined with a hollow cube (a 50 cube minus a 40 cube, so the 20 cube
  sits in the cavity), then cut by a 10 cube. The issue does not place the 10 cube.
- **Reported (SolidCraft):** the result fails SmBrep check #12, "Inner shell is not
  geometrically inside its containing outerShell". Both operands are valid.
- **Rust:** with the 10 cube overlapping the island's corner (10..20 against the island's
  15..35), the result is a manifold solid with V = 50³ - 40³ + 20³ - 5³. Other placements of
  the cutter are untested.
- **Test:** `crates/solidcore/tests/issue_388_island_in_void.rs`,
  `issue_388_island_in_a_void_stays_valid`.

These replace the draft upstream-issue records WR-362, WR-363 and WR-364, which were never
merged; those IDs are unused.
