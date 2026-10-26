//! Topology at a vertex all of whose edges are filleted (the NxN corners).
//!
//! Upstream: `source/SMLib/src/SmFilletCorner.cpp` (5196–5338, 5642–5852).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`FilletCorner::make_nxn_corner_topology`] | `MakeCornerTopology` of `SmFilletConvexNxNCorner` / `SmFilletConcaveNxNCorner` (by kind) |
//! | `make_convex_nxn_corner_topology` (private) | `SmFilletConvexNxNCorner::MakeCornerTopology` |
//! | `make_concave_nxn_corner_topology` (private) | `SmFilletConcaveNxNCorner::MakeCornerTopology` |
//! | [`NxNTopologyEnv`] | the corner's executive (`GetFilletSolverOfEdgeuse`) and pseudo-B-rep |
//!
//! Deviations: the caller lends the target, the pseudo-B-rep and the executive's
//! `GetFilletSolverOfEdgeuse` as a lookup (`None` for an edge use no solver
//! fillets). A missing solver or common face (`NER` upstream) is
//! `SM_ERR_NULL_POINTER`; topology upstream dereferences without a check is an input
//! error. The convex corner's unused `m_cpVertex->GetEdges` is not called. There is
//! no recursion.

use crate::blend::fillet_corner::p1::face_by_two_edgeuses;
use crate::blend::fillet_corner::p2::find_concave_edgeuses;
use crate::blend::fillet_corner::p3::{FilletCorner, FilletCornerKind};
use crate::blend::fillet_geom::p1::{
    CornerLookup, FilletBrep, FilletCornerId, FilletEdgeType, FilletGeomId, FilletVertex,
    FilletVertexType,
};
use crate::blend::fillet_solver::p1::FilletSolver;
use crate::error::{EntityRef, Result, SolidError, SM_ERR};
use crate::intersect::point_class::PointClassType;
use crate::kernel::Kernel;
use crate::math::Point3;
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::topology::edgeuse::p1 as eu1;
use crate::topology::{Brep, EdgeuseId, FaceId, Vertex, VertexId};

const OP: &str = "blend::fillet_corner::nxn_topology::make_nxn_corner_topology";

/// Actual objects used by the NxN corners' `MakeCornerTopology`.
pub struct NxNTopologyEnv<'a> {
    /// The executive's kernel.
    pub k: &'a Kernel,
    /// The target whose edge uses the corner stores.
    pub target: &'a mut Brep,
    /// The executive's actual pseudo-B-rep.
    pub pseudo: &'a mut FilletBrep,
    /// This corner's identity in the executive's corner store.
    pub corner: FilletCornerId,
    /// `GetFilletSolverOfEdgeuse`: the loaded solver filleting an edge use, `None` when
    /// none does.
    pub solver_of: &'a dyn Fn(EdgeuseId) -> Option<&'a FilletSolver>,
}

impl core::fmt::Debug for NxNTopologyEnv<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("NxNTopologyEnv")
            .field("corner", &self.corner)
            .finish_non_exhaustive()
    }
}

fn invalid(reason: &'static str) -> SolidError {
    SolidError::invalid_input(OP, reason)
}

fn sm_err() -> SolidError {
    SolidError::Upstream {
        operation: OP,
        code: SM_ERR,
        entity: None,
    }
}

fn need<T>(value: Option<T>, reason: &'static str) -> Result<T> {
    value.ok_or_else(|| invalid(reason))
}

/// `NER(p)`.
fn ner<T>(value: Option<T>) -> Result<T> {
    value.ok_or_else(|| SolidError::null_pointer(OP))
}

/// `GetFilletSolverOfEdgeuse(pEU)`, else of its mate (`NER` of the second).
fn solver_or_mate<'s>(
    target: &Brep,
    solver_of: &dyn Fn(EdgeuseId) -> Option<&'s FilletSolver>,
    eu: EdgeuseId,
) -> Result<&'s FilletSolver> {
    match solver_of(eu) {
        Some(solver) => Ok(solver),
        None => ner(solver_of(need(
            eu1::mate(target, eu)?,
            "edge use without mate",
        )?)),
    }
}

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
    /// `pFilletGeom->MakeRailEdgeuse(m_pPseudoBrep, lRailIndex, pFV, pNewEdgeuse)` with
    /// this corner as the running corner.
    #[allow(clippy::too_many_arguments)]
    fn nxn_rail_edgeuse(
        &self,
        pseudo: &mut FilletBrep,
        target: &Brep,
        corner: FilletCornerId,
        geometry: FilletGeomId,
        solver: &FilletSolver,
        rail: usize,
        fv: VertexId,
    ) -> Result<()> {
        pseudo.make_rail_edgeuse_on(
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
        Ok(())
    }

    /// `MakeCornerTopology` of an NxN corner (every edge at the vertex filleted): the
    /// convex corner (`SmFilletConvexNxNCorner`) puts a fillet vertex where each pair
    /// of neighbouring rails meets in their common face (two with a setback, joined by
    /// a rail-rail interpolation edge) and joins them with cross sections, or, for a
    /// three-edge bevel, with fillet/fillet edges to a centre vertex; the concave corner
    /// (`SmFilletConcaveNxNCorner`, exactly two concave edges) does the same with an
    /// extra vertex and a rail-rail interpolation edge in the face between its two
    /// concave edges, and defers setbacks to the convex corner. Sets `m_bBlending`.
    ///
    /// Allocations and corner-list insertions occur in upstream order and remain
    /// visible if a later operation fails.
    ///
    /// # Errors
    /// `SM_ERR` for a bevel of other than three edges or a concave corner without
    /// exactly two concave edges; `SM_ERR_NULL_POINTER` for a missing solver or common
    /// face; invalid topology or a corner of another kind; those of the convexity test,
    /// the solvers' rail lookups and the pseudo-topology operations.
    pub fn make_nxn_corner_topology(&mut self, env: NxNTopologyEnv<'_>) -> Result<()> {
        match self.kind {
            FilletCornerKind::ConvexNByN => self.make_convex_nxn_corner_topology(env),
            FilletCornerKind::ConcaveNByN => self.make_concave_nxn_corner_topology(env),
            _ => Err(invalid("requires an NxN corner")),
        }
    }

    /// `SmFilletConvexNxNCorner::MakeCornerTopology`.
    #[allow(clippy::too_many_lines)] // one upstream routine
    fn make_convex_nxn_corner_topology(&mut self, env: NxNTopologyEnv<'_>) -> Result<()> {
        let NxNTopologyEnv {
            k,
            target,
            pseudo,
            corner,
            solver_of,
        } = env;

        // Determine if this is a setback corner
        let is_set_back = self.set_back_dist > self.approx_tol_3d * 100.0;

        // All the edges of corner are filleted
        let all: Vec<EdgeuseId> = need(self.all_eus.clone(), "corner without outward edge uses")?;
        let total_eus = all.len();
        // (*m_pAllEUs)[lTotalEUs-1]
        let mut prev_eu = *need(all.last(), "corner without outward edge uses")?;
        let mut prev_solver = ner(solver_of(prev_eu))?;
        let last_solver = prev_solver;

        // Create all fillet vertices
        for &curr_eu in &all {
            let curr_solver = solver_or_mate(target, solver_of, curr_eu)?;
            let common_face: FaceId =
                ner(face_by_two_edgeuses(target, Some(prev_eu), Some(curr_eu), None)?)?;
            let prev_rail = prev_solver.find_index_of_rail_on_face(target, Some(common_face))?;
            let curr_rail = curr_solver.find_index_of_rail_on_face(target, Some(common_face))?;
            let curr_fv = new_corner_vertex(k, pseudo, corner)?;
            pseudo
                .fillet_vertex_mut(curr_fv)?
                .set_point_class(PointClassType::Face, Some(EntityRef::Face(common_face)));
            let prev_geometry = prev_solver.first_fillet_geom()?;
            let curr_geometry = curr_solver.first_fillet_geom()?;
            if is_set_back {
                // Will insert a new vertex before curr vert
                let new_v = new_corner_vertex(k, pseudo, corner)?;
                pseudo
                    .fillet_vertex_mut(new_v)?
                    .set_point_class(PointClassType::Face, Some(EntityRef::Face(common_face)));
                pseudo
                    .fillet_vertex_mut(new_v)?
                    .set_fillet_vertex_type(FilletVertexType::Mate);
                self.vertices.push(new_v);
                self.vertices.push(curr_fv);
                self.nxn_rail_edgeuse(
                    pseudo,
                    target,
                    corner,
                    prev_geometry,
                    prev_solver,
                    prev_rail,
                    new_v,
                )?;
                pseudo
                    .fillet_vertex_mut(curr_fv)?
                    .set_fillet_vertex_type(FilletVertexType::Setback);
                self.nxn_rail_edgeuse(
                    pseudo,
                    target,
                    corner,
                    curr_geometry,
                    curr_solver,
                    curr_rail,
                    curr_fv,
                )?;
            } else {
                // Regular corner, no setbacks
                pseudo
                    .fillet_vertex_mut(curr_fv)?
                    .set_fillet_vertex_type(FilletVertexType::RailXRail);
                self.nxn_rail_edgeuse(
                    pseudo,
                    target,
                    corner,
                    prev_geometry,
                    prev_solver,
                    prev_rail,
                    curr_fv,
                )?;
                self.nxn_rail_edgeuse(
                    pseudo,
                    target,
                    corner,
                    curr_geometry,
                    curr_solver,
                    curr_rail,
                    curr_fv,
                )?;
                self.vertices.push(curr_fv);
            }
            prev_solver = curr_solver;
            prev_eu = curr_eu;
        }

        // Define center vert for Bevel cases
        let mut center_fv = None;
        if self.bevel {
            self.blending = false;
            if total_eus != 3 {
                return Err(sm_err()); // Will not handle other cases for now
            }
            let center = new_corner_vertex(k, pseudo, corner)?;
            pseudo
                .fillet_vertex_mut(center)?
                .set_fillet_vertex_type(FilletVertexType::FilletX2Fillets);
            center_fv = Some(center);
        } else {
            // This is a blending corner (it blends all edge-fillets)
            self.blending = true;
        }

        // Now, create all fillet edges
        let total = self.vertices.len();
        prev_solver = last_solver;
        let mut count = 0_usize;
        let mut curr_fv = *need(self.vertices.first(), "corner without vertices")?;
        for &curr_eu in &all {
            let curr_solver = solver_or_mate(target, solver_of, curr_eu)?;
            let prev_geometry = prev_solver.first_fillet_geom()?;
            let curr_geometry = curr_solver.first_fillet_geom()?;
            // Get next vertex
            count = (count + 1) % total;
            let mut next_fv = *need(self.vertices.get(count), "corner vertex index")?;
            if is_set_back {
                // Insert a rail-rail-interpo edge
                let e = pseudo.make_fillet_edge(Some(curr_fv), Some(next_fv), Some(corner))?;
                pseudo
                    .fillet_edge_mut(e)?
                    .set_fillet_edge_type(FilletEdgeType::RailRailInterpolation);
                let common_face = match pseudo.fillet_vertex(next_fv)?.point_class_object() {
                    Some(EntityRef::Face(f)) => f,
                    _ => return Err(SolidError::null_pointer(OP)),
                };
                pseudo.fillet_edge_mut(e)?.set_original_face(Some(common_face));
                self.edges.push(e);
                curr_fv = next_fv;
                count = (count + 1) % total;
                next_fv = *need(self.vertices.get(count), "corner vertex index")?;
            }
            let fillet_edge = if let Some(center) = center_fv {
                // Make edge from curr vert to center
                let e = pseudo.make_fillet_edge(Some(curr_fv), Some(center), Some(corner))?;
                pseudo
                    .fillet_edge_mut(e)?
                    .set_fillet_edge_type(FilletEdgeType::FilletXFillet);
                let primary = need(
                    eu1::primary_edgeuse(pseudo.brep(), e)?,
                    "fillet edge without primary edge use",
                )?;
                let mate = need(
                    eu1::mate(pseudo.brep(), primary)?,
                    "fillet edge use without mate",
                )?;
                pseudo.set_edgeuse_fillet_geom(primary, Some(prev_geometry))?;
                pseudo.set_edgeuse_fillet_geom(mate, Some(curr_geometry))?;
                e
            } else {
                // Make a cross-section (iso-curve) on fillet
                let e = pseudo.make_fillet_edge(Some(curr_fv), Some(next_fv), Some(corner))?;
                pseudo
                    .fillet_edge_mut(e)?
                    .set_fillet_edge_type(FilletEdgeType::CrossSection);
                let primary = need(
                    eu1::primary_edgeuse(pseudo.brep(), e)?,
                    "fillet edge without primary edge use",
                )?;
                pseudo.set_edgeuse_fillet_geom(primary, Some(curr_geometry))?;
                if is_set_back {
                    pseudo.fillet_vertex_mut(curr_fv)?.set_mate(0, Some(next_fv));
                    pseudo.fillet_vertex_mut(next_fv)?.set_mate(0, Some(curr_fv));
                }
                e
            };
            self.edges.push(fillet_edge);
            curr_fv = next_fv;
            prev_solver = curr_solver;
        }
        // Add center vertex for BEVEL case
        if let Some(center) = center_fv {
            self.vertices.push(center);
        }
        Ok(())
    }

    /// `SmFilletConcaveNxNCorner::MakeCornerTopology`.
    #[allow(clippy::too_many_lines)] // one upstream routine
    fn make_concave_nxn_corner_topology(&mut self, env: NxNTopologyEnv<'_>) -> Result<()> {
        // Determine if this is a setback corner
        if self.set_back_dist > self.approx_tol_3d * 100.0 {
            // Let SmFilletConvexNxNCorner handle setbacks
            return self.make_convex_nxn_corner_topology(env);
        }
        let NxNTopologyEnv {
            k,
            target,
            pseudo,
            corner,
            solver_of,
        } = env;
        let vertex = need(self.vertex, "corner without vertex")?;
        let all: Vec<EdgeuseId> = need(self.all_eus.clone(), "corner without outward edge uses")?;

        // All the edges of corner are filleted
        let concave_eus = find_concave_edgeuses(k, target, vertex, &all)?;
        let &[concave_1, concave_2] = concave_eus.as_slice() else {
            return Err(sm_err());
        };
        let total_eus = all.len();
        let mut total = total_eus;

        // Allocate 'lTotalEUs' Fillet Vertices
        for _ in 0..total {
            let fv = new_corner_vertex(k, pseudo, corner)?;
            self.vertices.push(fv);
        }
        let mut prev_eu = *need(all.last(), "corner without outward edge uses")?;
        let mut prev_solver = ner(solver_of(prev_eu))?;
        let mut prev_fv = *need(self.vertices.last(), "corner vertex index")?;
        let mut jj = 0_usize;
        for &curr_eu in &all {
            let curr_fv = *need(self.vertices.get(jj), "corner vertex index")?;
            let common_face: FaceId =
                ner(face_by_two_edgeuses(target, Some(prev_eu), Some(curr_eu), None)?)?;
            pseudo
                .fillet_vertex_mut(curr_fv)?
                .set_point_class(PointClassType::Face, Some(EntityRef::Face(common_face)));
            let curr_solver = ner(solver_of(curr_eu))?;
            let rail = prev_solver.find_index_of_rail_on_face(target, Some(common_face))?;
            let other_rail = curr_solver.find_index_of_rail_on_face(target, Some(common_face))?;
            let prev_geometry = prev_solver.first_fillet_geom()?;
            let curr_geometry = curr_solver.first_fillet_geom()?;
            if (prev_eu == concave_1 && curr_eu == concave_2)
                || (prev_eu == concave_2 && curr_eu == concave_1)
            {
                // Process concave corner of pCommonFace: add one more fillet point
                let extra = new_corner_vertex(k, pseudo, corner)?;
                self.vertices.insert(jj, extra);
                jj += 1;
                total += 1;
                pseudo
                    .fillet_vertex_mut(extra)?
                    .set_point_class(PointClassType::Face, Some(EntityRef::Face(common_face)));
                pseudo
                    .fillet_vertex_mut(extra)?
                    .set_fillet_vertex_type(FilletVertexType::Mate);
                self.nxn_rail_edgeuse(
                    pseudo,
                    target,
                    corner,
                    prev_geometry,
                    prev_solver,
                    rail,
                    extra,
                )?;
                pseudo.fillet_vertex_mut(extra)?.set_mate(0, Some(prev_fv));
                pseudo.fillet_vertex_mut(prev_fv)?.set_mate(0, Some(extra));
                // Make fillet edge
                let e = pseudo.make_fillet_edge(Some(prev_fv), Some(extra), Some(corner))?;
                pseudo
                    .fillet_edge_mut(e)?
                    .set_fillet_edge_type(FilletEdgeType::CrossSection);
                let primary = need(
                    eu1::primary_edgeuse(pseudo.brep(), e)?,
                    "fillet edge without primary edge use",
                )?;
                pseudo.set_edgeuse_fillet_geom(primary, Some(prev_geometry))?;
                self.edges.push(e);
                // Setup curr vert
                pseudo
                    .fillet_vertex_mut(curr_fv)?
                    .set_fillet_vertex_type(FilletVertexType::Mate);
                self.nxn_rail_edgeuse(
                    pseudo,
                    target,
                    corner,
                    curr_geometry,
                    curr_solver,
                    other_rail,
                    curr_fv,
                )?;
                let next_fv = *need(
                    self.vertices.get((jj + 1) % total),
                    "corner vertex index",
                )?;
                pseudo.fillet_vertex_mut(curr_fv)?.set_mate(0, Some(next_fv));
                pseudo.fillet_vertex_mut(next_fv)?.set_mate(0, Some(curr_fv));
                // Make fillet edge between extra vert and current vert
                let e = pseudo.make_fillet_edge(Some(extra), Some(curr_fv), Some(corner))?;
                pseudo
                    .fillet_edge_mut(e)?
                    .set_fillet_edge_type(FilletEdgeType::RailRailInterpolation);
                pseudo.fillet_edge_mut(e)?.set_original_face(Some(common_face));
                self.edges.push(e);
            } else {
                pseudo
                    .fillet_vertex_mut(curr_fv)?
                    .set_fillet_vertex_type(FilletVertexType::RailXRail);
                // Setup rails of fillet
                self.nxn_rail_edgeuse(
                    pseudo,
                    target,
                    corner,
                    prev_geometry,
                    prev_solver,
                    rail,
                    curr_fv,
                )?;
                self.nxn_rail_edgeuse(
                    pseudo,
                    target,
                    corner,
                    curr_geometry,
                    curr_solver,
                    other_rail,
                    curr_fv,
                )?;
                // Make fillet edge
                let e = pseudo.make_fillet_edge(Some(prev_fv), Some(curr_fv), Some(corner))?;
                pseudo
                    .fillet_edge_mut(e)?
                    .set_fillet_edge_type(FilletEdgeType::CrossSection);
                let primary = need(
                    eu1::primary_edgeuse(pseudo.brep(), e)?,
                    "fillet edge without primary edge use",
                )?;
                pseudo.set_edgeuse_fillet_geom(primary, Some(prev_geometry))?;
                self.edges.push(e);
            }
            prev_solver = curr_solver;
            prev_fv = curr_fv;
            prev_eu = curr_eu;
            jj += 1;
        }
        self.blending = true;
        Ok(())
    }
}
