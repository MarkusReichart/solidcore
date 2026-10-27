## general: filleting a `Brep::clone()` fails in the merge (`MergeCurveOnSurfaces` 1006)

**What happens:** `fillet_executive::api::fillet_edges` on a `create_box` body works. On `body.clone()` the same call fails with `boolean::topology_intersector::p3::merge_curve_on_surfaces: upstream status 1006` (`SM_ERR_NULL_POINTER`). SolidCraft clones bodies for undo, history and previews. The SolidCraft cases runner hit it through its `take()`, which clones; its fillet step now moves the body out of the map instead.

**Cause:** `Brep::clone()` is upstream's copy constructor (`Brep::copy(false)`). By main's design it leaves every copied surface and curve without an owner, because `Surface`/`Curve` clones reset `owner`, as their C++ copy constructors do. Main tests that this is deliberate:
- `topology::tests::copy_standing_for_the_brep_keeps_geometry_owners` ("the copy constructor resets the owners");
- `tests/constructed_transforms.rs::clone_does_not_rebind_unowned_geometry` (0be3e8377: a transform of a clone leaves its unowned geometry in place).

The fillet merge resolves each merged surface to its face through the surface's owner (`InsertionGeometry::surface_face` → `owner_binding`). An owner-less surface resolves to no face, and `MergeCurveOnSurfaces` returns `SM_ERR_NULL_POINTER` (`topology_intersector/p3.rs`, the two `surface_face(..).ok_or(null_pointer)` sites).

**Upstream:** a body is copied with `SmBrep(const SmBrep&)` / `SmApiBrepCopy`, which re-owns every copied surface and curve with `SetOwner`. In the port that is `topology::brep::p3::copy_brep` / `api_brep_copy`, not `Clone`.

**Tried and dropped (fixer-7, 2026-10-26):** re-owning the geometry in `Clone` (restoring each owner, plus recording the source's names in `geometry_owner_keys`). It fixed the fillet but broke both main tests above, so `Clone` was left as main defines it.

**Open (for the owner to decide):** either
- (a) document that a body to be operated on must be copied with `api_brep_copy` (as upstream copies), not `clone()`, and give SolidCraft a public copy entry; or
- (b) let `clone()` keep owned geometry owned while leaving unowned geometry unowned. That would contradict `copy_standing_for_the_brep_keeps_geometry_owners`' plain-clone assertion, which checks that owned geometry's owner is reset.

Reproduce: `create_box` (100×100×10), `let mut c = b.clone();`, `fillet_edges(&mut k, &mut c, &[first edge], 1.0, 1, 1, 1.0)` gives 1006.
