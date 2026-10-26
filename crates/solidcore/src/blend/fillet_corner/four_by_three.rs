//! Topology of the 4x3 fillet corner: four edges at the corner, three filleted.
//!
//! Upstream: `source/SMLib/src/SmFilletCorner.cpp` (443–506, 8570–8669).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`FilletCorner::make_4x3_corner_topology`] | `SmFillet4x3Corner::MakeCornerTopology` |
//! | [`FourByThreeTopologyEnv`] | the executive's pseudo-B-rep and `GetFilletSolverOfEdgeuse` |
//! | `find_side_edgeuse` (private) | `sm_FindSideEdgeuse` (as `topology.rs`'s) |
//!
//! Deviations: the caller lends the target, the pseudo-B-rep and the executive's
//! solver lookup ([`FourByThreeTopologyEnv`]). Missing topology that upstream dereferences is
//! an input error; a corner without filleted edge uses (upstream divides by zero) too.
//! There is no recursion.

use crate::blend::fillet_corner::p1::face_by_two_edgeuses;
use crate::blend::fillet_corner::p3::{FilletCorner, FilletCornerKind};
use crate::blend::fillet_geom::p1::{
    CornerLookup, FilletBrep, FilletCornerId, FilletEdgeType, FilletGeomId, FilletVertex,
    FilletVertexType,
};
use crate::blend::fillet_solver::p1::FilletSolver;
use crate::error::{EntityRef, Result, SolidError};
use crate::intersect::point_class::PointClassType;
use crate::kernel::Kernel;
use crate::math::Point3;
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::topology::edgeuse::p1 as eu1;
use crate::topology::{Brep, EdgeuseId, Vertex, VertexId};

const OP: &str = "blend::fillet_corner::four_by_three::make_4x3_corner_topology";

/// What `SmFillet4x3Corner::MakeCornerTopology` works on.
pub struct FourByThreeTopologyEnv<'a> {
    /// The executive's kernel.
    pub k: &'a Kernel,
    /// The B-rep being filleted.
    pub target: &'a Brep,
    /// The executive's pseudo-B-rep.
    pub pseudo: &'a mut FilletBrep,
    /// This corner's identity in the executive's corner store.
    pub corner: FilletCornerId,
    /// `GetFilletSolverOfEdgeuse`: the loaded solver of an edge use, if any.
    pub solver_of: &'a dyn Fn(EdgeuseId) -> Option<&'a FilletSolver>,
}

impl core::fmt::Debug for FourByThreeTopologyEnv<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("FourByThreeTopologyEnv")
            .field("corner", &self.corner)
            .finish_non_exhaustive()
    }
}

fn invalid(reason: &'static str) -> SolidError {
    SolidError::invalid_input(OP, reason)
}

fn need<T>(value: Option<T>, reason: &'static str) -> Result<T> {
    value.ok_or_else(|| invalid(reason))
}

/// `NER(p)`.
fn ner<T>(value: Option<T>) -> Result<T> {
    value.ok_or_else(|| SolidError::null_pointer(OP))
}

/// `sm_FindSideEdgeuse(pFilEU, pSideEU)`: the side edge use or its radial, the one
/// whose loop neighbour is the filleted edge.
fn find_side_edgeuse(
    target: &Brep,
    filleted: EdgeuseId,
    side: EdgeuseId,
) -> Result<Option<EdgeuseId>> {
    let filleted_edge = eu1::edge(target, filleted)?;
    let cw = need(
        eu1::cw_edgeuse(target, side)?,
        "side edge use without CW neighbour",
    )?;
    if filleted_edge == eu1::edge(target, cw)? {
        return Ok(Some(side));
    }
    let side = need(eu1::radial(target, side)?, "edge use without radial")?;
    let ccw = need(
        eu1::ccw_edgeuse(target, side)?,
        "side radial without CCW neighbour",
    )?;
    Ok((filleted_edge == eu1::edge(target, ccw)?).then_some(side))
}

/// `new (m_pPseudoBrep) SmFilletVertex(this)`.
fn new_corner_vertex(
    k: &Kernel,
    pseudo: &mut FilletBrep,
    corner: FilletCornerId,
) -> Result<VertexId> {
    let vertex = pseudo.brep_mut().add_unlinked_vertex_record(Vertex {
        point: Point3::ZERO,
        tolerance: 0.0,
        ..Vertex::default()
    })?;
    pseudo.insert_fillet_vertex(vertex, FilletVertex::new(k, Some(corner)));
    Ok(vertex)
}

impl FilletCorner {
    /// `pFilletGeom->MakeRailEdgeuse(m_pPseudoBrep, rail, fv, pNewEdgeuse)`.
    fn rail_4x3(
        &self,
        env: &mut FourByThreeTopologyEnv<'_>,
        geom: FilletGeomId,
        solver: &FilletSolver,
        rail: usize,
        fv: VertexId,
    ) -> Result<()> {
        env.pseudo.make_rail_edgeuse_on(
            geom,
            rail,
            Some(fv),
            OrientType::Unknown,
            None,
            CornerLookup {
                this: Some((Some(env.corner), self)),
                ..Default::default()
            },
            solver,
            env.target,
        )?;
        Ok(())
    }

    /// `SmFillet4x3Corner::MakeCornerTopology`: four edges at the corner, three of them
    /// filleted. Around the corner, each filleted edge use gets a vertex where its
    /// fillet's rail meets the next fillet's: in their common face, or (after the
    /// unfilleted edge use) on the side edge, where the previous and next fillets' rails
    /// cross it. Consecutive vertices are joined by cross sections of the next fillet.
    ///
    /// # Errors
    ///
    /// Invalid topology or a different corner kind; `NER` for a missing solver or
    /// common face; those of the pseudo-topology operations.
    pub fn make_4x3_corner_topology(&mut self, mut env: FourByThreeTopologyEnv<'_>) -> Result<()> {
        if self.kind != FilletCornerKind::FourByThree {
            return Err(invalid("requires a 4x3 corner"));
        }
        let all = need(self.all_eus.clone(), "corner without outward edge uses")?;
        let total_eus = all.len();
        let total_solver_eus = need(self.solver_eus.as_ref(), "corner without solver edge uses")?
            .len();
        if total_eus == 0 || total_solver_eus == 0 {
            return Err(invalid("corner without edge uses"));
        }

        // Allocate the fillet vertices
        for _ in 0..total_solver_eus {
            let fv = new_corner_vertex(env.k, env.pseudo, env.corner)?;
            self.vertices.push(fv);
        }

        let solver_of = env.solver_of;
        let mut curr_eu = *need(all.first(), "corner without edge uses")?;
        let mut curr_solver = solver_of(curr_eu);
        let mut curr_index = 0;
        let mut curr_fv = need(self.vertices.first().copied(), "fillet vertex index")?;
        for ii in 0..total_eus {
            let next_eu = *need(all.get((ii + 1) % total_eus), "edge use index")?;
            let Some(next_solver) = solver_of(next_eu) else {
                curr_solver = None;
                curr_eu = next_eu;
                continue;
            };
            let next_geom = next_solver.first_fillet_geom()?;

            if let Some(curr) = curr_solver {
                let common_face = ner(face_by_two_edgeuses(
                    env.target,
                    Some(curr_eu),
                    Some(next_eu),
                    None,
                )?)?;
                let fv = env.pseudo.fillet_vertex_mut(curr_fv)?;
                fv.set_point_class(PointClassType::Face, Some(EntityRef::Face(common_face)));
                fv.set_fillet_vertex_type(FilletVertexType::RailXRail);
                // Create rail edgeuses
                let curr_geom = curr.first_fillet_geom()?;
                let rail = curr.find_index_of_rail_on_face(env.target, Some(common_face))?;
                let other = next_solver.find_index_of_rail_on_face(env.target, Some(common_face))?;
                self.rail_4x3(&mut env, curr_geom, curr, rail, curr_fv)?;
                self.rail_4x3(&mut env, next_geom, next_solver, other, curr_fv)?;
            } else {
                let prev_eu = *need(
                    all.get((ii + total_eus - 1) % total_eus),
                    "edge use index",
                )?;
                let prev_solver = ner(solver_of(prev_eu))?;
                let prev_geom = prev_solver.first_fillet_geom()?;
                let side = need(
                    find_side_edgeuse(env.target, prev_eu, curr_eu)?,
                    "no side edge use intersects the rail",
                )?;
                let fv = env.pseudo.fillet_vertex_mut(curr_fv)?;
                fv.set_point_class(PointClassType::Edgeuse, Some(EntityRef::Edgeuse(side)));
                fv.set_fillet_vertex_type(FilletVertexType::RailXEdgeuse);
                // Create rail edgeuses
                let rail = prev_solver.find_index_of_rail_x_side_edgeuse(env.target, side)?;
                let radial = need(eu1::radial(env.target, side)?, "edge use without radial")?;
                let other = next_solver.find_index_of_rail_x_side_edgeuse(env.target, radial)?;
                self.rail_4x3(&mut env, prev_geom, prev_solver, rail, curr_fv)?;
                self.rail_4x3(&mut env, next_geom, next_solver, other, curr_fv)?;
            }

            // Make fillet edge
            curr_index += 1;
            let next_fv = need(
                self.vertices.get(curr_index % total_solver_eus).copied(),
                "fillet vertex index",
            )?;
            let edge = env.pseudo.make_fillet_edge(Some(curr_fv), Some(next_fv), Some(env.corner))?;
            env.pseudo
                .fillet_edge_mut(edge)?
                .set_fillet_edge_type(FilletEdgeType::CrossSection);
            let primary = need(
                eu1::primary_edgeuse(env.pseudo.brep(), edge)?,
                "fillet edge without primary edge use",
            )?;
            env.pseudo.set_edgeuse_fillet_geom(primary, Some(next_geom))?;
            self.edges.push(edge);
            curr_solver = Some(next_solver);
            curr_eu = next_eu;
            curr_fv = next_fv;
        }
        Ok(())
    }
}
