//! Rail classification and rollover processing for a fillet solver.
//!
//! Upstream: `source/SMLib/src/SmFilletSolver.cpp` (5714–5822, 7161–7899).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`FilletSolver::test_roll_over`] | `SmFilletSolver::TestRollOver` |
//! | [`FilletSolver::process_roll_over`] | `SmFilletSolver::ProcessRollOver` |
//! | [`FilletSolver::rails_fit_their_faces`] | `SmFilletSolver::RailsFitTheirFaces` (SolidCraft fork) |
//! | `leaves_past_an_end` (private) | `SmFilletSolver::LeavesPastAnEnd` (SolidCraft fork) |
//!
//! Deviations: pseudo-B-rep objects use typed arenas and geometry uses the same
//! typed ids as the solver list. Completed mutations survive failures. Only unavailable
//! callees are stubbed; the entire caller continuation is ported, including both
//! full-rollover calls, the segment walk and topology/surface/corner updates.
//! The four-argument classification constructor calls the real `reset`
//! (`intersect.curve_class_p2`), preserving its defaults and source tolerance.
//! Internal short intervals are retained in descending order for refinement.
//! See `p5/recovery.rs` for the bounded segment-walk deviation. No recursion.
//! Opt-in (WR-301, `Kernel::fixes.fillet_rail_rim_on_face`): before the short-interval
//! squeeze, neighbouring intervals of a rail's classification on its face or on the
//! face's boundary are joined, so only the parts of a rail outside its face roll over
//! (upstream takes a rail running along its face's rim as leaving the face and fails
//! "Unable to Process Rollover"). Off by default: the C++ comparison is unchanged.
//! Not upstream (WR-302): [`FilletSolver::rails_fit_their_faces`] ports the SolidCraft
//! fork's rail check (tangency-labs/smlib f0edce8ce, with `LeavesPastAnEnd` from
//! tangency-labs/solidcore#841 and the off-face run of #1024). It only answers a
//! question; the executive calls it behind the opt-in
//! `Kernel::fixes.fillet_refuse_rails_leaving_faces`.

use super::p1::FilletSolver;
use crate::blend::fillet_classification::{FilletClassObject, LegacyClassificationObjects};
use crate::blend::fillet_geom::p1::{
    FilletBrep, FilletCornerId, FilletCornerType, FilletGeomId, FilletGeomType,
};
#[cfg(test)]
use crate::curve::curve::p2::OwnerChange;
use crate::error::{ignored, EntityRef, Result, SolidError};
#[cfg(test)]
use crate::geom::object::ObjectRef;
use crate::intersect::curve_class::p1::{CurveClassification, IntervalPosition};
use crate::intersect::curve_class::p2::ResetOptions;
use crate::intersect::point_class::{ClassObject, PointClassType};
use crate::kernel::Kernel;
use crate::query::face::p2::{curve_on_classify, CurveOnClassifyObjects};
use crate::surface::surf_of_revolution::SurfOfRevolution;
use crate::topology::{edge, edgeuse::p1 as eu, edgeuse::p2, Brep, EdgeId, EdgeuseId};

#[path = "p5/recovery.rs"]
mod recovery;

const OP: &str = "blend::fillet_solver::p5";

/// Resolve handles used by `TestRollOver`, without assuming the two B-reps' arena
/// indices identify the same objects. Curves returned here must be resolvable by
/// the inherited classification geometry resolver.
pub trait RolloverObjects: CurveOnClassifyObjects<EntityRef> {
    /// Identity of a rail's 3D curve in the pseudo-B-rep.
    fn rail_curve(&self, rail: EdgeId) -> EntityRef;
    /// Identity of the original face's UV curve on the rail edgeuse's mate.
    fn rail_mate_uv_curve(&self, mate: EdgeuseId) -> Option<EntityRef>;
    /// The actual virtual `GetCornerType` result for the specified corner.
    ///
    /// # Errors
    /// A resolver must report an unresolved corner, rather than invent its type.
    fn corner_type(&self, corner: FilletCornerId) -> Result<FilletCornerType>;
}

/// Actual identity and classification services, borrowed only during each classification.
pub trait RailClassificationContext {
    /// Actual retained rail curve identity.
    fn rail_curve(&self, rail: EdgeId) -> FilletClassObject;
    /// Actual retained mate UV identity, when present.
    fn rail_mate_uv_curve(&self, pseudo: &Brep, mate: EdgeuseId) -> Option<FilletClassObject>;
    /// Classify the attached record using the two actual topology domains.
    ///
    /// # Errors
    /// Classifier errors, retaining all partial record writes.
    fn classify(
        &mut self,
        k: &Kernel,
        target: &mut Brep,
        pseudo: &Brep,
        face: crate::topology::FaceId,
        class: &mut CurveClassification<FilletClassObject>,
    ) -> Result<()>;
    /// Actual corner kind.
    ///
    /// # Errors
    /// An unresolved corner.
    fn corner_type(&self, corner: FilletCornerId) -> Result<FilletCornerType>;
}

struct LegacyRolloverContext<'a, T>(&'a T);
impl<T: RolloverObjects> RailClassificationContext for LegacyRolloverContext<'_, T> {
    fn rail_curve(&self, rail: EdgeId) -> FilletClassObject {
        FilletClassObject::Legacy(self.0.rail_curve(rail))
    }
    fn rail_mate_uv_curve(&self, _: &Brep, mate: EdgeuseId) -> Option<FilletClassObject> {
        self.0
            .rail_mate_uv_curve(mate)
            .map(FilletClassObject::Legacy)
    }
    fn classify(
        &mut self,
        k: &Kernel,
        target: &mut Brep,
        _: &Brep,
        face: crate::topology::FaceId,
        class: &mut CurveClassification<FilletClassObject>,
    ) -> Result<()> {
        curve_on_classify(
            k,
            target,
            face,
            &LegacyClassificationObjects(self.0),
            true,
            class,
            true,
            false,
            true,
        )?;
        Ok(())
    }
    fn corner_type(&self, corner: FilletCornerId) -> Result<FilletCornerType> {
        self.0.corner_type(corner)
    }
}

pub use super::p1::RailOwnership as RolloverOwnership;

impl FilletSolver {
    /// Classify both rails against their original faces, squeeze tiny end
    /// intervals, refine their classifications and process rollover. Iterates a
    /// snapshot of the solver's geometry list, as upstream does. Every entry is
    /// the geometry's actual id in the caller's pseudo-B-rep arena.
    ///
    /// # Errors
    /// Missing rails/curves/mates are invalid inputs (upstream dereferences them).
    /// Classification and unported recovery errors propagate; completed earlier
    /// changes and the caller's adjusted-corner list are retained.
    pub fn test_roll_over(
        &mut self,
        k: &Kernel,
        target: &mut Brep,
        fillets: &mut FilletBrep,
        objects: &impl RolloverObjects,
        adjusted_corners: &mut Vec<FilletCornerId>,
        ownership: &mut impl RolloverOwnership,
    ) -> Result<()> {
        self.test_roll_over_with_context(
            k,
            target,
            fillets,
            &mut LegacyRolloverContext(objects),
            adjusted_corners,
            ownership,
        )
    }

    /// Classify actual retained rails through a scoped, short-lived geometry provider.
    ///
    /// # Errors
    /// As [`Self::test_roll_over`], preserving attached partial classifications.
    pub fn test_roll_over_with_context(
        &mut self,
        k: &Kernel,
        target: &mut Brep,
        fillets: &mut FilletBrep,
        context: &mut impl RailClassificationContext,
        adjusted_corners: &mut Vec<FilletCornerId>,
        ownership: &mut impl RolloverOwnership,
    ) -> Result<()> {
        let geoms = self.fillet_geoms();
        for id in geoms {
            for ii in 0..2 {
                let rail = fillets
                    .geom(id)?
                    .rail(ii)
                    .ok_or_else(|| SolidError::invalid_input(OP, "TestRollOver: missing rail"))?;
                let curve = edge::curve(fillets.brep(), rail)?.ok_or_else(|| {
                    SolidError::invalid_input(OP, "TestRollOver: missing rail curve")
                })?;
                let interval = curve.natural_interval();
                let primary = eu::primary_edgeuse(fillets.brep(), rail)?.ok_or_else(|| {
                    SolidError::invalid_input(OP, "TestRollOver: missing primary edgeuse")
                })?;
                let mate = eu::mate(fillets.brep(), primary)?
                    .ok_or_else(|| SolidError::invalid_input(OP, "TestRollOver: missing mate"))?;
                let mut cc = CurveClassification::new(k);
                cc.reset(
                    k,
                    &ResetOptions {
                        curve: Some(context.rail_curve(rail)),
                        interval: Some(interval),
                        uv_curve: context.rail_mate_uv_curve(fillets.brep(), mate),
                        zone_tol_3d: self.approx_tol_3d,
                        ..ResetOptions::default()
                    },
                )?;
                fillets.fillet_edge_mut(rail)?.set_curve_class(Some(cc));
                let face = fillets.fillet_edge(rail)?.original_face().ok_or_else(|| {
                    SolidError::invalid_input(OP, "TestRollOver: missing original face")
                })?;
                // Install the classification before calling the fallible operation:
                // upstream leaves that object attached even when classification fails.
                fillets.with_rail_classification(rail, |pseudo, cc| {
                    context.classify(k, target, pseudo, face, cc)
                })?;
                finish_classification(
                    self,
                    k,
                    target,
                    fillets,
                    id,
                    ii,
                    rail,
                    &mut MissingClassification,
                )?;
            }
            self.process_roll_over(
                k,
                target,
                fillets,
                id,
                &|corner| context.corner_type(corner),
                adjusted_corners,
                ownership,
            )?;
        }
        Ok(())
    }

    /// Process an existing pair of rail classifications. Blend patches and
    /// single-interval rails that both remain classified need no repair.
    /// `corner_type` resolves the corner on rail 0's start vertex only when full
    /// rollover is needed; `adjusted_corners` is append-only upstream output.
    ///
    /// # Errors
    /// A missing classification or corner returns `SM_ERR_NULL_POINTER`. Unsupported
    /// full-rollover corner types return `SM_ERR`. Recovery through unavailable
    /// units returns `Unported`; invalid handles are reported by the arena.
    pub fn process_roll_over(
        &mut self,
        k: &Kernel,
        target: &mut Brep,
        fillets: &mut FilletBrep,
        id: FilletGeomId,
        corner_type: &impl Fn(FilletCornerId) -> Result<FilletCornerType>,
        adjusted_corners: &mut Vec<FilletCornerId>,
        ownership: &mut impl RolloverOwnership,
    ) -> Result<()> {
        recovery::Recovery {
            solver: self,
            k,
            target,
            fillets,
            ownership,
        }
        .process(id, corner_type, adjusted_corners, &mut recovery::Missing)
    }
}

impl FilletSolver {
    /// Not upstream (the SolidCraft fork's `SmFilletSolver::RailsFitTheirFaces`, WR-302):
    /// whether every rail of this solver's default fillet geometries stays on the face
    /// it rolls on, so that the blend fits. Each untrimmed rail is classified against
    /// its original face (`SmFace::CurveOnClassify`, as `TestRollOver`). It may run off
    /// the face past the filleted edge's ends (`LeavesPastAnEnd`; the corner trims it
    /// there) and off across a tangent side edge (a rollover onto the next face). A rail
    /// nowhere on its face, or leaving it beside the edge across an edge that is not
    /// tangent (a wall, an opening's rim), means the blend is larger than the faces
    /// allow: `false`. A run of off-face intervals (split where the rail passes through
    /// a boundary point from outside to outside, tangency-labs/solidcore#1024) leaves
    /// the face once and is judged at its first interval.
    ///
    /// # Errors
    /// A missing rail, rail curve, original face, edge use or mate (the fork's `NER`),
    /// and those of the classification, the tangency test and `LeavesPastAnEnd`.
    pub fn rails_fit_their_faces(
        &self,
        k: &Kernel,
        target: &mut Brep,
        fillets: &FilletBrep,
        objects: &impl RolloverObjects,
    ) -> Result<bool> {
        self.rails_fit_their_faces_with_context(
            k,
            target,
            fillets,
            &mut LegacyRolloverContext(objects),
        )
    }

    /// [`Self::rails_fit_their_faces`] through a scoped classification provider.
    ///
    /// # Errors
    /// As [`Self::rails_fit_their_faces`].
    pub fn rails_fit_their_faces_with_context(
        &self,
        k: &Kernel,
        target: &mut Brep,
        fillets: &FilletBrep,
        context: &mut impl RailClassificationContext,
    ) -> Result<bool> {
        let null = || SolidError::null_pointer(OP);
        for id in self.fillet_geoms() {
            let geom = fillets.geom(id)?;
            if geom.geom_type != FilletGeomType::Default {
                continue;
            }
            for ii in 0..2 {
                let rail = geom.rail(ii).ok_or_else(null)?;
                let curve = edge::curve(fillets.brep(), rail)?.ok_or_else(null)?.clone();
                let face = fillets.fillet_edge(rail)?.original_face().ok_or_else(null)?;
                let primary = eu::primary_edgeuse(fillets.brep(), rail)?.ok_or_else(null)?;
                let mate = eu::mate(fillets.brep(), primary)?.ok_or_else(null)?;
                let mut class = CurveClassification::new(k);
                class.reset(
                    k,
                    &ResetOptions {
                        curve: Some(context.rail_curve(rail)),
                        interval: Some(curve.natural_interval()),
                        uv_curve: context.rail_mate_uv_curve(fillets.brep(), mate),
                        zone_tol_3d: self.approx_tol_3d,
                        ..ResetOptions::default()
                    },
                )?;
                context.classify(k, target, fillets.brep(), face, &mut class)?;

                let (on_face, runs) = off_face_runs(&class)?;
                let mut fit = true;
                for run in runs {
                    // The rail beyond the face from there to its end: past an end of the
                    // filleted edge, where a corner trims it, only if it leaves the face
                    // there.
                    if run.to_rail_end
                        && self.leaves_past_an_end(k, target, face, &curve, run.param, &run.leave)?
                    {
                        continue;
                    }
                    let side = match target_object(&run.leave, PointClassType::Edge) {
                        Some(EntityRef::Edge(e)) => Some(e),
                        _ => None,
                    };
                    let tangent = match side {
                        Some(e) => {
                            super::p1::test_manifold_tangent_rollover(
                                k,
                                target,
                                e,
                                self.tangency_tolerance,
                            )?
                            .0
                        }
                        None => false,
                    };
                    if !tangent {
                        fit = false;
                    }
                }
                if !on_face || !fit {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    /// Not upstream (the SolidCraft fork's `SmFilletSolver::LeavesPastAnEnd`,
    /// tangency-labs/solidcore#841): whether a rail leaves `face` at `leave` (rail
    /// parameter `rail_param`) past an end of the filleted edge, where the corner there
    /// trims it. Either across the face's side edge at an end (the edge before or after
    /// the filleted edge in the face's loop) or through one of that edge's vertices; or
    /// on or beyond the cross section at an end: the filleted edge's point nearest to
    /// where it leaves is an end, within the solver's approximation tolerance.
    fn leaves_past_an_end(
        &self,
        k: &Kernel,
        target: &Brep,
        face: crate::topology::FaceId,
        rail: &crate::curve::curve::Curve,
        rail_param: f64,
        leave: &crate::intersect::point_class::PointClassification<FilletClassObject>,
    ) -> Result<bool> {
        use crate::geom::global_solver::SolverOperation;
        let null = || SolidError::null_pointer(OP);
        let fil_eu = self.edgeuse(0)?.ok_or_else(null)?;
        let fil_edge = eu::edge(target, fil_eu)?.ok_or_else(null)?;
        let edge_curve = edge::curve(target, fil_edge)?.ok_or_else(null)?;

        // Across or through the face's side edge at an end.
        let leave_edge = match target_object(leave, PointClassType::Edge) {
            Some(EntityRef::Edge(e)) => Some(e),
            _ => None,
        };
        let leave_vertex = match target_object(leave, PointClassType::Vertex) {
            Some(EntityRef::Vertex(v)) => Some(v),
            _ => None,
        };
        for edge_use in edge::edgeuses(target, fil_edge, None, None)? {
            if eu::face(target, edge_use)? != Some(face)
                || !eu::is_loop_edgeuse(target, edge_use)?
            {
                continue;
            }
            let cw = eu::cw_edgeuse(target, edge_use)?;
            let ccw = eu::ccw_edgeuse(target, edge_use)?;
            for next in [cw, ccw] {
                let side = match next {
                    Some(next) => eu::edge(target, next)?,
                    None => None,
                };
                let Some(side) = side.filter(|&s| s != fil_edge) else {
                    continue;
                };
                let v0 = edge::vertex(target, side)?;
                let v1 = edge::other_vertex(target, side, v0)?;
                if leave_edge == Some(side)
                    || (leave_vertex.is_some() && (leave_vertex == v0 || leave_vertex == v1))
                {
                    return Ok(true);
                }
            }
        }

        // On or beyond the cross section at an end.
        let edge_ivl = edge::interval(target, fil_edge)?;
        let leave_point = rail.evaluate_point(k, rail_param)?;
        let dropped = edge_curve.drop_point(
            k,
            &edge_ivl,
            leave_point,
            None,
            crate::math::BIG_DOUBLE,
            None,
            SolverOperation::Minimize,
        )?;
        if !dropped.success {
            return Err(SolidError::sm_err(OP));
        }
        let foot = edge_curve.evaluate_point(k, dropped.parameter)?;
        for end in [edge_ivl.min(), edge_ivl.max()] {
            let end = edge_curve.evaluate_point(k, end)?;
            if foot.distance(end) <= self.approx_tol_3d {
                return Ok(true);
            }
        }
        Ok(false)
    }
}

/// A run of a rail's classification intervals off its face (`SM_PC_UNKNOWN` mids), for
/// `RailsFitTheirFaces`.
#[derive(Debug, Clone, Copy)]
struct OffFaceRun<O> {
    /// Where the rail leaves the face: the end of the run when it starts the rail, else
    /// the start of its first interval.
    leave: crate::intersect::point_class::PointClassification<O>,
    /// The rail parameter there.
    param: f64,
    /// Whether the run reaches the rail's start or end (the rail beyond the face from
    /// `leave` to that end).
    to_rail_end: bool,
}

/// The fork's interval walk of `RailsFitTheirFaces`: whether some interval of `class`
/// is on the face (or its boundary), and each run of intervals off it, judged once at
/// its first interval (tangency-labs/solidcore#1024: a run split where the rail passes
/// through a boundary point from outside to outside leaves the face once).
fn off_face_runs<O: ClassObject>(
    class: &CurveClassification<O>,
) -> Result<(bool, Vec<OffFaceRun<O>>)> {
    let off = |j: usize| -> Result<bool> {
        Ok(class.at(j)?.mid.point_class() == PointClassType::Unknown)
    };
    let num_ivls = class.size();
    let mut on_face = false;
    let mut runs = Vec::new();
    for j in 0..num_ivls {
        let ivl = *class.at(j)?;
        if !off(j)? {
            on_face = true;
            continue;
        }
        if j > 0 && off(j - 1)? {
            continue;
        }
        let mut j_last = j;
        while j_last + 1 < num_ivls && off(j_last + 1)? {
            j_last += 1;
        }
        // Where the rail leaves the face: the first interval ends there, the others start.
        let last = *class.at(j_last)?;
        let (leave, param) = if j == 0 {
            (last.end, last.interval.max())
        } else {
            (ivl.start, ivl.interval.min())
        };
        runs.push(OffFaceRun {
            leave,
            param,
            to_rail_end: j == 0 || j_last + 1 == num_ivls,
        });
    }
    Ok((on_face, runs))
}

/// The target object of a point classification of class `class` (`GetObject()` cast
/// to the class's type), `None` for another class or a non-target object.
fn target_object(
    pc: &crate::intersect::point_class::PointClassification<FilletClassObject>,
    class: PointClassType,
) -> Option<EntityRef> {
    if pc.point_class() != class {
        return None;
    }
    pc.object().and_then(|o| o.target_entity().ok())
}

trait ClassificationOps {
    fn connect(
        &mut self,
        _k: &Kernel,
        fillets: &mut FilletBrep,
        rail: EdgeId,
        interval: usize,
        at_end: bool,
        allow_removal_of_ends: bool,
    ) -> Result<()> {
        fillets
            .fillet_edge_mut(rail)?
            .curve_class
            .as_mut()
            .ok_or_else(|| SolidError::null_pointer(OP))?
            .connect_intervals(
                interval,
                if at_end {
                    IntervalPosition::End
                } else {
                    IntervalPosition::Start
                },
                allow_removal_of_ends,
            )
    }
    fn refine(
        &mut self,
        solver: &mut FilletSolver,
        k: &Kernel,
        target: &mut Brep,
        fillets: &mut FilletBrep,
        _geom: FilletGeomId,
        which: usize,
        rail: EdgeId,
        small: &[usize],
    ) -> Result<()> {
        let mut classification = fillets
            .fillet_edge_mut(rail)?
            .curve_class
            .take()
            .ok_or_else(|| SolidError::null_pointer(OP))?;
        let result = super::p3::refine_rail_curve_classification_scoped(
            solver,
            k,
            target,
            fillets,
            which,
            small,
            rail,
            &mut classification,
        );
        // The classification is detached only to share immutable rail data with
        // the callee; restore its partial writes on every return path.
        fillets.fillet_edge_mut(rail)?.curve_class = Some(classification);
        result
    }
}
struct MissingClassification;
impl ClassificationOps for MissingClassification {}

fn finish_classification(
    solver: &mut FilletSolver,
    k: &Kernel,
    target: &mut Brep,
    fillets: &mut FilletBrep,
    geom: FilletGeomId,
    which: usize,
    rail: EdgeId,
    ops: &mut impl ClassificationOps,
) -> Result<()> {
    if k.fixes.fillet_rail_rim_on_face {
        let class = fillets
            .fillet_edge_mut(rail)?
            .curve_class
            .as_mut()
            .ok_or_else(|| SolidError::null_pointer(OP))?;
        join_on_face_intervals(class)?;
    }
    let small = squeeze_rail(k, fillets, rail, solver.approx_tol_3d, ops)?;
    // Refinement is warning-only upstream, but retain its completed mutations
    // and propagate an unavailable implementation rather than faking success.
    ignored(ops.refine(solver, k, target, fillets, geom, which, rail, &small))?;
    Ok(())
}

/// Not upstream (WR-301, opt-in `Kernel::fixes.fillet_rail_rim_on_face`): the rail leaves
/// its face only where it runs outside it. Where it touches the face's boundary at a
/// vertex or runs along a piece of an edge of it (a blend whose setback reaches
/// exactly to the rim of a hole), it is still on the face: neighbouring intervals
/// whose mids are on the face or its boundary become one interval
/// (`ConnectIntervals(jj, SM_IP_START, FALSE)`, from the last pair down), taking the
/// later mid where the earlier one is not `SM_PC_FACE`.
fn join_on_face_intervals<O: ClassObject>(class: &mut CurveClassification<O>) -> Result<()> {
    let mut jj = class.size();
    while jj > 1 {
        jj -= 1;
        let prev_mid = class.at(jj - 1)?.mid;
        let this_mid = class.at(jj)?.mid;
        if prev_mid.point_class() == PointClassType::Unknown
            || this_mid.point_class() == PointClassType::Unknown
        {
            continue;
        }
        if prev_mid.point_class() != PointClassType::Face {
            class.at_mut(jj - 1)?.mid = this_mid;
        }
        class.connect_intervals(jj, IntervalPosition::Start, false)?;
    }
    Ok(())
}

fn squeeze_rail(
    k: &Kernel,
    fillets: &mut FilletBrep,
    rail: EdgeId,
    tolerance: f64,
    ops: &mut impl ClassificationOps,
) -> Result<Vec<usize>> {
    let count = fillets
        .fillet_edge(rail)?
        .curve_class()
        .ok_or_else(|| SolidError::null_pointer(OP))?
        .size();
    let mut small = Vec::new();
    if count > 1 {
        for index in (0..count).rev() {
            let interval = fillets
                .fillet_edge(rail)?
                .curve_class()
                .ok_or_else(|| SolidError::null_pointer(OP))?
                .sub_interval(index)?;
            let curve =
                edge::curve(fillets.brep(), rail)?.ok_or_else(|| SolidError::null_pointer(OP))?;
            let start = curve.evaluate_point(k, interval.min())?;
            let end = curve.evaluate_point(k, interval.max())?;
            if start.distance(end) > 50.0 * tolerance {
                continue;
            }
            if index == 0 {
                ops.connect(k, fillets, rail, 0, false, true)?;
            } else if index == count - 1 {
                ops.connect(k, fillets, rail, index, true, true)?;
            } else {
                small.push(index);
            }
        }
    }
    Ok(small)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blend::fillet_geom::p1::{FilletEdge, FilletGeom, FilletVertex};
    use crate::intersect::point_class::PointClassType;
    use crate::math::Vec3;

    /// tangency-labs/solidcore#485: a rail classified on the face, on its boundary edge
    /// and on the face again (FACE / EDGE / FACE, the ends on vertices) is one interval
    /// on the face with the opt-in (WR-301); intervals outside the face stay apart.
    #[test]
    fn issue_485_rail_along_rim_stays_on_face() {
        let k = Kernel::default();
        let face = Some(EntityRef::Face(crate::topology::FaceId::from_index(0)));
        let edge = Some(EntityRef::Edge(EdgeId::from_index(3)));
        let classes = [
            PointClassType::Face,
            PointClassType::Edge,
            PointClassType::Face,
            PointClassType::Unknown,
            PointClassType::Edge,
            PointClassType::Face,
        ];
        let mut class: CurveClassification<EntityRef> = CurveClassification::new(&k);
        class.set_size(&k, classes.len()).unwrap();
        for (i, &pc) in classes.iter().enumerate() {
            let ivl = class.at_mut(i).unwrap();
            ivl.interval.set_min_max(i as f64, (i + 1) as f64).unwrap();
            ivl.start.set_class_object(PointClassType::Vertex, None);
            ivl.end.set_class_object(PointClassType::Vertex, None);
            let object = if pc == PointClassType::Edge { edge } else { face };
            ivl.mid.set_class_object(pc, object);
        }
        let mut upstream = class.clone();
        // upstream: every boundary stays an interval boundary
        assert_eq!(upstream.size(), 6);
        join_on_face_intervals(&mut class).unwrap();
        let mids: Vec<_> = class
            .intervals()
            .iter()
            .map(|ivl| (ivl.interval.min(), ivl.interval.max(), ivl.mid.point_class()))
            .collect();
        assert_eq!(
            mids,
            [
                (0.0, 3.0, PointClassType::Face),
                (3.0, 4.0, PointClassType::Unknown),
                (4.0, 6.0, PointClassType::Face),
            ]
        );
        // a single interval is left as it is
        upstream.set_size(&k, 1).unwrap();
        join_on_face_intervals(&mut upstream).unwrap();
        assert_eq!(upstream.size(), 1);
    }

    /// A classification with one interval per `(mid, start, end)` class, the interval
    /// `[i, i + 1]`.
    fn classified(
        k: &Kernel,
        classes: &[(PointClassType, PointClassType, PointClassType)],
    ) -> CurveClassification<EntityRef> {
        let edge = Some(EntityRef::Edge(EdgeId::from_index(3)));
        let mut class: CurveClassification<EntityRef> = CurveClassification::new(k);
        class.set_size(k, classes.len()).unwrap();
        for (i, &(mid, start, end)) in classes.iter().enumerate() {
            let ivl = class.at_mut(i).unwrap();
            ivl.interval.set_min_max(i as f64, (i + 1) as f64).unwrap();
            // (no object makes a class Unknown)
            let object = if mid == PointClassType::Unknown { None } else { edge };
            ivl.mid.set_class_object(mid, object);
            ivl.start.set_class_object(start, edge);
            ivl.end.set_class_object(end, edge);
        }
        class
    }

    /// tangency-labs/solidcore#1024 (in the fork's `RailsFitTheirFaces`, WR-302): a rail
    /// that runs off its face through a boundary vertex has two adjacent off-face
    /// intervals there; the run leaves the face once, where its first interval starts.
    /// A run from the rail's start leaves where its last interval ends.
    #[test]
    fn issue_1024_off_face_run_through_a_vertex_is_judged_once() {
        use PointClassType::{Edge, Face, Unknown, Vertex};
        let k = Kernel::default();
        let class = classified(
            &k,
            &[
                (Face, Vertex, Edge),
                (Unknown, Edge, Vertex),
                (Unknown, Vertex, Edge),
                (Face, Edge, Vertex),
            ],
        );
        let (on_face, runs) = off_face_runs(&class).unwrap();
        assert!(on_face);
        assert_eq!(runs.len(), 1);
        assert_eq!(runs[0].param, 1.0);
        assert_eq!(runs[0].leave.point_class(), Edge);
        assert!(!runs[0].to_rail_end);

        let class = classified(
            &k,
            &[
                (Unknown, Vertex, Vertex),
                (Unknown, Vertex, Edge),
                (Face, Edge, Vertex),
                (Unknown, Edge, Vertex),
            ],
        );
        let (on_face, runs) = off_face_runs(&class).unwrap();
        assert!(on_face);
        let leaves: Vec<_> = runs
            .iter()
            .map(|r| (r.param, r.leave.point_class(), r.to_rail_end))
            .collect();
        assert_eq!(leaves, [(2.0, Edge, true), (3.0, Edge, true)]);

        let class = classified(&k, &[(Unknown, Vertex, Vertex)]);
        let (on_face, runs) = off_face_runs(&class).unwrap();
        assert!(!on_face);
        assert_eq!(runs.len(), 1);
    }

    struct NoOwnership;
    impl RolloverOwnership for NoOwnership {
        fn edge_owner(&self, _: EdgeId) -> Result<(ObjectRef, Option<ObjectRef>)> {
            Err(SolidError::invalid_input(OP, "unexpected owner lookup"))
        }
        fn owner_changed(&mut self, _: OwnerChange) -> Result<()> {
            Err(SolidError::invalid_input(
                OP,
                "unexpected owner notification",
            ))
        }
    }

    fn classified_rails(k: &Kernel, fb: &mut FilletBrep) -> FilletGeomId {
        let mut geom = FilletGeom::new(k, None, FilletGeomType::Default).unwrap();
        for i in 0..2 {
            let a = fb.brep_mut().add_vertex(Vec3::ZERO).unwrap();
            let b = fb.brep_mut().add_vertex(Vec3::X).unwrap();
            let edge = fb.brep_mut().add_edge(a, b).unwrap();
            let vu = fb
                .brep_mut()
                .add_vertexuse_record(crate::topology::Vertexuse {
                    vertex: Some(a),
                    ..crate::topology::Vertexuse::default()
                })
                .unwrap();
            let eu = fb
                .brep_mut()
                .add_edgeuse_record(crate::topology::Edgeuse {
                    edge: Some(edge),
                    vertexuse: Some(vu),
                    ..crate::topology::Edgeuse::default()
                })
                .unwrap();
            fb.brep_mut().edge_mut(edge).unwrap().edgeuses.push(eu);

            fb.insert_fillet_vertex(a, FilletVertex::new(k, Some(FilletCornerId(42 + i as u32))));
            let mut record = FilletEdge::with_corner(None);
            let mut cc = CurveClassification::new(k);
            cc.set_size(k, 1).unwrap();
            cc.at_mut(0)
                .unwrap()
                .mid
                .set_class_object(PointClassType::Edge, Some(EntityRef::Edge(edge)));
            record.set_indexed_curve_class(Some(&cc)).unwrap();
            fb.insert_fillet_edge(edge, record);
            geom.set_rail(i, Some(edge));
        }
        fb.add_geom(geom)
    }

    #[test]
    fn resolver_is_lazy_and_failure_preserves_caller_output() {
        let k = Kernel::default();
        let mut fb = FilletBrep::new();
        let id = classified_rails(&k, &mut fb);
        let mut solver = FilletSolver::new(1e-5, 1e-3, 1e-3);
        let mut target = Brep::new();
        let mut adjusted = vec![FilletCornerId(9), FilletCornerId(9)];
        solver
            .process_roll_over(
                &k,
                &mut target,
                &mut fb,
                id,
                &|_| panic!("known rails must not resolve a corner"),
                &mut adjusted,
                &mut NoOwnership,
            )
            .unwrap();
        let rail = fb.geom(id).unwrap().rail(1).unwrap();
        fb.fillet_edge_mut(rail)
            .unwrap()
            .curve_class
            .as_mut()
            .unwrap()
            .at_mut(0)
            .unwrap()
            .mid
            .set_class_object(PointClassType::Unknown, None);
        let refusal = SolidError::invalid_input(OP, "unresolved corner");
        let result = solver.process_roll_over(
            &k,
            &mut target,
            &mut fb,
            id,
            &|corner| {
                assert_eq!(corner, FilletCornerId(42));
                Err(refusal.clone())
            },
            &mut adjusted,
            &mut NoOwnership,
        );
        assert_eq!(result, Err(refusal));
        assert_eq!(adjusted, [FilletCornerId(9), FilletCornerId(9)]);
        let result = solver.process_roll_over(
            &k,
            &mut target,
            &mut fb,
            id,
            &|_| Ok(FilletCornerType::OneByOne),
            &mut adjusted,
            &mut NoOwnership,
        );
        assert_eq!(
            result,
            Err(SolidError::invalid_input(
                "blend::fillet_solver::p3",
                "solver without edge use",
            ))
        );
    }
}
