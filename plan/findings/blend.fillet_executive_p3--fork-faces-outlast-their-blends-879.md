# blend.fillet_executive_p3: the SolidCraft fork's `FacesOutlastTheirBlends` (#879) is an opt-in fork feature, not an upstream defect (2026-10-26)

Unit: `blend.fillet_executive_p3` (`SmFilletExecutive::CreateFilletBrep`). PR: the swarm batch
carrying tangency-labs/solidcore#879.

Upstream SMLib has no `FacesOutlastTheirBlends`: the check is a function of the SolidCraft
fork only, so a blend that takes the whole of a face it rolls on is upstream's designed
behaviour, not an NVIDIA defect. This record replaces the former
`plan/upstream-issues/WR-320--smfilletexecutive-createfilletbrep-builds-blends-that-take-the.md`
(withdrawn: fork-only, per the 2026-10-26 upstream-issues audit). The port keeps the check as
an opt-in fork feature; default behaviour is SMLib's.

- **Kind:** fork feature (behaviour SolidCraft's fork adds; SMLib builds the blend)
- **Area:** Blends (fillets, chamfers)
- **Location:** `source/SMLib/src/SmFilletExecutive.cpp` `SmFilletExecutive::CreateFilletBrep`: nothing between the rails' trimming and classification and the fillet faces' construction asks whether a face the blends roll on keeps any of itself.
- **What happens:** a blend as wide as the face it rolls on is built. A radius-10, height-20 cylinder's top rim chamfered 10 across the top comes back as a 3-face body: the top face is gone and the chamfer's cone ends in a point at the top's centre. The rail on the top has shrunk to that point (every point of its B-spline at (0, 0, 20)), and a point rail by itself is fine (the corners of tangency-labs/solidcore#420 have them), so nothing notices that no part of the top is left.
- **Trigger:** every edge of a loop of a face blended, each of those blends' rails on that face a point, and some blend there with a cross section that meets the face at an angle (the linear section of a chamfer). A round's arc meets the face tangentially at its rail, so a round as wide as the top closes the end as a hemisphere, tangent to the side and smooth at the pole; that is a valid result.
- **How found:** tangency-labs/solidcore#879 (the SolidCraft fork's `smlib_chamfer_sides` and its kernel tests), with gdb on the fork.
- **Confirmed:** in the fork's C++ build (the 3-face body passes SMLib's level-2 check); not run here.
- **Rust port:** Faithful by default: `FilletExecutive::create_fillet_brep` builds it as upstream does. With the opt-in `Kernel::fixes.fillet_refuse_face_consuming_blends` on, the check of the fork's `SmFilletExecutive::FacesOutlastTheirBlends` (with its amendment for tangent sections) runs after `classify_solver_rails`: a loop of a face that closes to a point, some blend on it not tangent to the face (`FilletSurfaceGenerator::meets_its_faces_tangentially`: base and custom `false`, linear `false`, circular `!mirror`, blend curve `true`), notes `FilletErrorType::NoRoom` ("the blends take the whole of a face they roll on, some at an angle to it") and fails with `SM_ERR`. The live session's `SessionObjects` answers the check (`BrepCreationAccess::creation_faces_outlast_their_blends`); adapters without the rails' geometry keep the default, no refusal. Tests: `crates/solidcore/tests/issue_879_face_consuming_blends.rs` (`issue_879_chamfer_taking_the_whole_top_is_refused`, `issue_879_round_taking_the_whole_top_passes_the_rule`, `issue_879_a_chamfer_short_of_the_whole_top_passes_the_rule`, `issue_879_default_keeps_upstream`).

Refs tangency-labs/solidcore#879
