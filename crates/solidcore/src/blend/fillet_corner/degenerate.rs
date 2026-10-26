//! Topology of the degenerate fillet corner: the fillet's rails meet at the corner's
//! own vertex.
//!
//! Upstream: `source/SMLib/src/SmFilletCorner.cpp` (8779–8861).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`FilletCorner::make_degenerate_corner_topology`] | `SmFilletDegenerateCorner::MakeCornerTopology` |
//! | [`DegenerateTopologyEnv`] | the executive's pseudo-B-rep and `GetFilletSolverOfEdgeuse` |
//!
//! Deviations: the caller lends the target, the pseudo-B-rep and the loaded solver of
//! the filleted edge use. The corner point drops onto each base surface with
//! `GlobalPointSolve` over its natural domain without point testing: upstream turns
//! the point test of the surface's cache off around the call; the port solves on an
//! owner-less evaluation copy, whose solve classifies nothing. Missing topology that
//! upstream dereferences is an input error. There is no recursion.

use crate::blend::fillet_corner::p3::{FilletCorner, FilletCornerKind};
use crate::blend::fillet_geom::p1::{
    CornerLookup, FilStatus, FilletBrep, FilletCornerId, FilletVertex, FilletVertexType,
};
use crate::blend::fillet_solver::p1::FilletSolver;
use crate::error::{EntityRef, Result, SolidError};
use crate::geom::global_solver::{SolutionRequested, SolverOperation};
use crate::intersect::point_class::PointClassType;
use crate::intersect::surface_intersector::p1::{IntersectionPointType, SrfIndex};
use crate::kernel::Kernel;
use crate::math::{Point2, Point3};
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::surface::surface::p4::NoObject;
use crate::topology::edgeuse::p1 as eu1;
use crate::topology::{Brep, Vertex};

const OP: &str = "blend::fillet_corner::degenerate::make_degenerate_corner_topology";

/// What `SmFilletDegenerateCorner::MakeCornerTopology` works on.
#[derive(Debug)]
pub struct DegenerateTopologyEnv<'a> {
    /// The executive's kernel.
    pub k: &'a Kernel,
    /// The B-rep being filleted.
    pub target: &'a Brep,
    /// The executive's pseudo-B-rep.
    pub pseudo: &'a mut FilletBrep,
    /// This corner's identity in the executive's corner store.
    pub corner: FilletCornerId,
    /// The loaded solver of the filleted edge use (`GetFilletSolverOfEdgeuse`).
    pub solver: &'a FilletSolver,
}

fn invalid(reason: &'static str) -> SolidError {
    SolidError::invalid_input(OP, reason)
}

fn need<T>(value: Option<T>, reason: &'static str) -> Result<T> {
    value.ok_or_else(|| invalid(reason))
}

impl FilletCorner {
    /// `SmFilletDegenerateCorner::MakeCornerTopology`: one fillet vertex on the corner's
    /// own vertex (at its point, already processed), where the rails crossing the
    /// filleted edge use's clockwise neighbour and its radial's counter-clockwise
    /// neighbour both end. The corner point drops onto each of the solver's base
    /// surfaces; each rail end records the two surfaces' parameters there as a tangent
    /// point.
    ///
    /// # Errors
    ///
    /// Invalid topology or a different corner kind; `NER` for a solver without its
    /// offset surfaces; `SM_ERR` unless the corner point drops onto each base surface
    /// at exactly one point; those of the point solve and the pseudo-topology
    /// operations.
    pub fn make_degenerate_corner_topology(
        &mut self,
        env: DegenerateTopologyEnv<'_>,
    ) -> Result<()> {
        let DegenerateTopologyEnv {
            k,
            target,
            pseudo,
            corner,
            solver,
        } = env;
        if self.kind != FilletCornerKind::Degenerate {
            return Err(invalid("requires a degenerate corner"));
        }
        let filleted = need(
            self.solver_eus
                .as_deref()
                .and_then(|eus| eus.first())
                .copied(),
            "corner without filleted edge use",
        )?;
        // the adjacent edge uses of the filleted edge use
        let radial = need(eu1::radial(target, filleted)?, "edge use without radial")?;
        let adjacent = [
            need(eu1::cw_edgeuse(target, filleted)?, "edge use without CW neighbour")?,
            need(eu1::ccw_edgeuse(target, radial)?, "radial without CCW neighbour")?,
        ];
        let vertex = need(self.vertex, "corner without vertex")?;
        let vertex_record = target.vertex(vertex)?;
        let (corner_pnt, vertex_tol) = (vertex_record.point, vertex_record.tolerance);

        let fv = pseudo.brep_mut().add_unlinked_vertex_record(Vertex {
            point: Point3::ZERO,
            tolerance: 0.0,
            ..Vertex::default()
        })?;
        pseudo.insert_fillet_vertex(fv, FilletVertex::new(k, Some(corner)));
        let record = pseudo.fillet_vertex_mut(fv)?;
        record.set_point_class(PointClassType::Vertex, Some(EntityRef::Vertex(vertex)));
        record.set_fillet_vertex_type(FilletVertexType::OnVertex);
        pseudo.set_vertex_point(fv, corner_pnt)?;
        pseudo.fillet_vertex_mut(fv)?.set_status(FilStatus::Processed);
        self.vertices.push(fv);

        // Drop the corner point onto the 2 surfaces
        let mut uv = [Point2::new(0.0, 0.0); 2];
        for (ii, slot) in uv.iter_mut().enumerate() {
            let offset = solver
                .surface(ii)?
                .ok_or_else(|| SolidError::null_pointer(OP))?;
            let base = offset.surface.base_surface();
            // point testing off: the solve keeps every point solution without
            // classifying it against the trim boundaries
            let mut standalone = base.clone_for_evaluation()?;
            let _ = standalone.set_owner(None);
            let solutions = standalone.global_point_solve::<NoObject>(
                k,
                Some(NoObject),
                &base.natural_uv_domain(),
                SolverOperation::Intersect,
                corner_pnt,
                vertex_tol * 10.0,
                None,
                SolutionRequested::Single,
            )?;
            let [solution] = solutions.as_slice() else {
                return Err(SolidError::sm_err(OP));
            };
            *slot = Point2::new(solution.start.param(0)?, solution.start.param(1)?);
        }

        let geometry = solver.first_fillet_geom()?;
        for eu in adjacent {
            let rail = solver.find_index_of_rail_x_side_edgeuse(target, eu)?;
            // Create rail-edgeuse using this vertex
            let new_eu = pseudo.make_rail_edgeuse_on(
                geometry,
                rail,
                Some(fv),
                OrientType::Unknown,
                None,
                CornerLookup {
                    this: Some((Some(corner), self)),
                    ..Default::default()
                },
                solver,
                target,
            )?;
            let vu = need(
                eu1::vertexuse(pseudo.brep(), new_eu)?,
                "rail edge use without vertex use",
            )?;
            let tsect = pseudo.fillet_vertexuse_mut(vu)?.tsect_pnt_mut();
            for side in [rail, 1 - rail] {
                let srf = if side == 0 {
                    SrfIndex::First
                } else {
                    SrfIndex::Second
                };
                *tsect.uv_pos_mut(srf) = need(uv.get(side).copied(), "rail index")?;
            }
            tsect.point_type = IntersectionPointType::TangentPoint;
        }
        Ok(())
    }
}
