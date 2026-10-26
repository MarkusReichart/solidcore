//! Topology at a vertex where two filleted edges meet (the Nx2 corner).
//!
//! Upstream: `source/SMLib/src/SmFilletCorner.cpp` (443–506, 548–604, 6487–6866,
//! 6905–7064).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`FilletCorner::make_nx2_corner_topology`] | `SmFilletNx2Corner::MakeCornerTopology` |
//! | [`Nx2TopologyEnv`] | the corner's executive (`GetFilletSolverOfEdgeuse`) and pseudo-B-rep |
//! | `make_4x2_concave_corner_topology` (private) | `SmFilletNx2Corner::Make4x2ConcaveCornerTopology` |
//! | `find_side_edgeuse` (private) | `sm_FindSideEdgeuse` (as `topology.rs`'s) |
//! | `find_angle_between_edgeuses` (private) | `sm_FindAngleBetweenEdgeuses` |
//!
//! Deviations: the caller lends the target, the pseudo-B-rep and the executive's
//! `GetFilletSolverOfEdgeuse` in an [`Nx2TopologyEnv`] (a missing solver of a
//! filleted edge use is `SM_ERR_NULL_POINTER`). Missing topology that upstream dereferences (no side
//! edge use found, no solver edge uses) is an input error. `sm_FindAngleBetweenEdgeuses`
//! returns a failed `GetEndTangent`'s status as the angle, as upstream's `SER` in a
//! `double` function does (`Unported` still propagates); a failed `AngleBetween`, whose
//! angle upstream leaves unset, is an input error. There is no recursion.

use crate::blend::fillet_corner::p1::{
    check_face_complimentary_corner_concavity, face_by_two_edgeuses, test_g1_continuity,
};
use crate::blend::fillet_corner::p3::{FilletCorner, FilletCornerKind};
use crate::blend::fillet_geom::p1::{
    CornerLookup, FilletBrep, FilletCornerId, FilletEdgeType, FilletGeomId, FilletVertex,
    FilletVertexType,
};
use crate::blend::fillet_solver::p1::FilletSolver;
use crate::error::{EntityRef, Result, SolidError, SM_ERR};
use crate::intersect::point_class::PointClassType;
use crate::kernel::Kernel;
use crate::math::{rad_to_deg, Point3};
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::topology::edge;
use crate::topology::edgeuse::p1 as eu1;
use crate::topology::{Brep, EdgeId, EdgeuseId, FaceId, Vertex, VertexId};

const OP: &str = "blend::fillet_corner::nx2_topology::make_nx2_corner_topology";

/// Actual objects used by `SmFilletNx2Corner::MakeCornerTopology`.
pub struct Nx2TopologyEnv<'a> {
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

impl core::fmt::Debug for Nx2TopologyEnv<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("Nx2TopologyEnv")
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

fn edge_of(target: &Brep, eu: EdgeuseId) -> Result<EdgeId> {
    need(eu1::edge(target, eu)?, "edge use without edge")
}

fn radial(target: &Brep, eu: EdgeuseId) -> Result<EdgeuseId> {
    need(eu1::radial(target, eu)?, "edge use without radial")
}

fn vertex_of(target: &Brep, eu: EdgeuseId) -> Result<Option<VertexId>> {
    let vu = need(eu1::vertexuse(target, eu)?, "edge use without vertex use")?;
    Ok(target.vertexuse(vu)?.vertex)
}

/// `sm_FindSideEdgeuse(pFilletedEU, pSideEU)`: the side edge use (or its radial) whose
/// neighbour around the corner is the filleted edge, `None` when neither is.
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
    let side = radial(target, side)?;
    let ccw = need(
        eu1::ccw_edgeuse(target, side)?,
        "side radial without CCW neighbour",
    )?;
    Ok((filleted_edge == eu1::edge(target, ccw)?).then_some(side))
}

/// Whether `vertex`'s end of `eu`'s edge is the edge's start (`bAtEdgeStart` of
/// `sm_FindAngleBetweenEdgeuses`).
fn at_edge_start(target: &Brep, vertex: VertexId, eu: EdgeuseId) -> Result<bool> {
    let at_vertex = vertex_of(target, eu)? == Some(vertex);
    let orientation = eu1::orientation(target, eu)?;
    Ok(!((at_vertex && orientation == OrientType::Opposite)
        || (!at_vertex && orientation == OrientType::Same)))
}

/// `GetEndTangent(bAtEdgeStart, sTangent)` reversed at the edge's end, so it points
/// away from the corner; `Err(angle)` is upstream's `SER` returning the status code as
/// `sm_FindAngleBetweenEdgeuses`'s `double`.
fn outward_tangent(
    k: &Kernel,
    target: &Brep,
    vertex: VertexId,
    eu: EdgeuseId,
) -> Result<core::result::Result<crate::math::Vec3, f64>> {
    let at_start = at_edge_start(target, vertex, eu)?;
    let e = edge_of(target, eu)?;
    match edge::end_tangent(k, target, e, at_start) {
        Ok(tangent) => Ok(Ok(if at_start { tangent } else { -tangent })),
        Err(e) if e.is_unported() => Err(e),
        Err(SolidError::Upstream { code, .. }) => Ok(Err(f64::from(code))),
        Err(e) => Err(e),
    }
}

/// `sm_FindAngleBetweenEdgeuses(cpCornerVert, pEU1, pEU2)`: the angle (radians)
/// between the two edges' tangents pointing away from the corner.
fn find_angle_between_edgeuses(
    k: &Kernel,
    target: &Brep,
    vertex: VertexId,
    eu_1: EdgeuseId,
    eu_2: EdgeuseId,
) -> Result<f64> {
    let tangent_1 = match outward_tangent(k, target, vertex, eu_1)? {
        Ok(t) => t,
        Err(status) => return Ok(status),
    };
    let tangent_2 = match outward_tangent(k, target, vertex, eu_2)? {
        Ok(t) => t,
        Err(status) => return Ok(status),
    };
    tangent_1
        .angle_between(tangent_2)
        .map_err(|_| invalid("AngleBetween failed (angle unset upstream)"))
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
    fn nx2_rail_edgeuse(
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

    /// `SmFilletNx2Corner::MakeCornerTopology`: the fillet vertices and edges of a
    /// corner where two of the three or four edges at the vertex are filleted. Each
    /// side (unfilleted) edge gets a fillet vertex where both rails cross it; with three
    /// edges a second vertex in the face the filleted edges share (the rails'
    /// intersection for a bevel; straight across from the side vertex for a tangent or
    /// shallow mixed-convexity corner; with a three-sided corner face in between for a
    /// mixed-convexity corner), with four the corner is a bevel or a blend between the
    /// two side vertices (a mixed-convexity side edge makes the 4x2 concave topology
    /// instead). Sets `m_bBevel`, `m_bBlending` and `m_bG1Case`.
    ///
    /// Allocations and corner-list insertions occur in upstream order and remain
    /// visible if a later operation fails.
    ///
    /// # Errors
    /// `SM_ERR` for more than four edges at the vertex, a filleted edge use missing from
    /// the corner's list, two adjacent filleted edge uses of four, or not exactly two
    /// vertices before the cross edge; invalid topology or a different corner kind;
    /// those of the tangency and convexity tests, the solvers' rail lookups and the
    /// pseudo-topology operations.
    #[allow(clippy::too_many_lines)] // one upstream routine
    pub fn make_nx2_corner_topology(&mut self, env: Nx2TopologyEnv<'_>) -> Result<()> {
        let Nx2TopologyEnv {
            k,
            target,
            pseudo,
            corner,
            solver_of,
        } = env;
        if self.kind != FilletCornerKind::NByTwo {
            return Err(invalid("requires an NByTwo corner"));
        }
        let all: Vec<EdgeuseId> = need(self.all_eus.clone(), "corner without outward edge uses")?;
        let total = all.len();

        // no work - unhandled case - Currently we handle only 3x2 and 4x2
        if total > 4 {
            return Err(sm_err());
        }

        // Nx2 corner, get the two filleted edgeuses from the base brep
        let solver_eus = need(self.solver_eus.as_deref(), "corner without solver edge uses")?;
        let base_1 = *need(solver_eus.first(), "corner without first filleted edge use")?;
        let base_2 = *need(solver_eus.get(1), "corner without second filleted edge use")?;
        let vertex = need(self.vertex, "corner without vertex")?;
        // GetFilletExecutive()->GetFilletSolverOfEdgeuse of the two filleted edge uses
        let solvers = [
            solver_of(base_1).ok_or_else(|| SolidError::null_pointer(OP))?,
            solver_of(base_2).ok_or_else(|| SolidError::null_pointer(OP))?,
        ];

        // for 4 EU case
        if total == 4 {
            // check state: to be filleted edgeuses are on list of outward edgeuses
            let (Some(index_1), Some(index_2)) = (
                all.iter().position(|&eu| eu == base_1),
                all.iter().position(|&eu| eu == base_2),
            ) else {
                return Err(sm_err());
            };

            // the two filleted edge uses must not be adjacent to each other
            if (index_1 + index_2) % 2 != 0 {
                return Err(sm_err());
            }

            // Check the non-filleted (side) edges for mixed convexity
            for &side in &all {
                if side == base_1 || side == base_2 {
                    continue; // (not a side edge)
                }
                if check_face_complimentary_corner_concavity(
                    k,
                    target,
                    vertex,
                    base_1,
                    base_2,
                    Some(side),
                )? {
                    return self.make_4x2_concave_corner_topology(
                        k, target, pseudo, corner, solvers, &all, base_1, base_2, vertex,
                    );
                }
            }
            self.bevel = true;
        } // end 4 outward EUs from this corner vertex case check

        // Determine if two edges are tangent at corner vertex.
        let tangency_tol_deg = rad_to_deg(self.tangency_tol_radians);
        let (is_tangent_case, angle_deg) = test_g1_continuity(
            k,
            target,
            Some(vertex),
            Some(edge_of(target, base_1)?),
            Some(edge_of(target, base_2)?),
            tangency_tol_deg,
        )?;
        // (upstream leaves dAngleDeg unset when the tangents cannot be compared; it is
        // read only below, where an unset angle is not a number)
        let angle_deg = angle_deg.unwrap_or(f64::NAN);

        // Store this value for later use.
        self.g1_case = is_tangent_case;

        let mixed_convexity =
            check_face_complimentary_corner_concavity(k, target, vertex, base_1, base_2, None)?;

        // Get Solvers and FilletGeoms for the two filleted edges at this vertex
        let [solver_1, solver_2] = solvers;
        let geometry_1 = solver_1.first_fillet_geom()?;
        let geometry_2 = solver_2.first_fillet_geom()?;

        // For the non-filleted edge(s), we expect the rails to hit it (them).
        // Create a (fillet) vertex there, along with two rail edgeuses.
        for &current in &all {
            if current == base_1 || current == base_2 {
                continue;
            }
            // This edgeuse is in between two filleted edgeuses
            let side = need(
                find_side_edgeuse(target, base_1, current)?,
                "no side edge use beside the first filleted edge",
            )?;
            let fv = new_corner_vertex(k, pseudo, corner)?;
            {
                let v = pseudo.fillet_vertex_mut(fv)?;
                v.set_point_class(PointClassType::Edgeuse, Some(EntityRef::Edgeuse(side)));
                v.set_fillet_vertex_type(FilletVertexType::RailXEdgeuse);
            }
            self.vertices.push(fv);

            let rail = solver_1.find_index_of_rail_x_side_edgeuse(target, side)?;
            let other_rail =
                solver_2.find_index_of_rail_x_side_edgeuse(target, radial(target, side)?)?;
            self.nx2_rail_edgeuse(pseudo, target, corner, geometry_1, solver_1, rail, fv)?;
            self.nx2_rail_edgeuse(pseudo, target, corner, geometry_2, solver_2, other_rail, fv)?;
        }

        // Create the topology for the fillet corner.
        if total == 3 {
            // Find the side EU -- the un-filleted one.
            let side = all.iter().copied().find(|&eu| eu != base_1 && eu != base_2);

            let common_face: Option<FaceId> =
                face_by_two_edgeuses(target, Some(base_1), Some(base_2), side)?;
            let rail = solver_1.find_index_of_rail_on_face(target, common_face)?;
            let other_rail = solver_2.find_index_of_rail_on_face(target, common_face)?;
            let side_fv = *need(self.vertices.first(), "corner without side vertex")?;
            let face_fv = new_corner_vertex(k, pseudo, corner)?;
            pseudo
                .fillet_vertex_mut(face_fv)?
                .set_point_class(PointClassType::Face, common_face.map(EntityRef::Face));
            self.vertices.push(face_fv);

            // mixed convexity case
            if mixed_convexity && !self.bevel && angle_deg > 5.0 {
                // Mixed-convexity 3x2 corner, will create a 3-sided corner face.
                self.blending = true;

                // straight across from where the other rail hits the side edge
                pseudo
                    .fillet_vertex_mut(face_fv)?
                    .set_fillet_vertex_type(FilletVertexType::Mate);
                pseudo.fillet_vertex_mut(face_fv)?.set_mate(0, Some(side_fv));
                pseudo.fillet_vertex_mut(side_fv)?.set_mate(0, Some(face_fv));

                self.nx2_rail_edgeuse(pseudo, target, corner, geometry_1, solver_1, rail, face_fv)?;

                // the vertex where the other rail meets the 3-sided face
                let face_fv_2 = new_corner_vertex(k, pseudo, corner)?;
                self.vertices.push(face_fv_2);
                pseudo
                    .fillet_vertex_mut(face_fv_2)?
                    .set_point_class(PointClassType::Face, common_face.map(EntityRef::Face));
                pseudo
                    .fillet_vertex_mut(face_fv_2)?
                    .set_fillet_vertex_type(FilletVertexType::Mate);
                pseudo.fillet_vertex_mut(face_fv_2)?.set_mate(0, Some(side_fv));
                pseudo.fillet_vertex_mut(side_fv)?.set_mate(1, Some(face_fv_2));

                self.nx2_rail_edgeuse(
                    pseudo, target, corner, geometry_2, solver_2, other_rail, face_fv_2,
                )?;

                // a straight-across rib to end filletGeom1
                let rib = pseudo.make_fillet_edge(Some(side_fv), Some(face_fv), Some(corner))?;
                pseudo
                    .fillet_edge_mut(rib)?
                    .set_fillet_edge_type(FilletEdgeType::CrossSection);
                let primary = need(
                    eu1::primary_edgeuse(pseudo.brep(), rib)?,
                    "rib without primary edge use",
                )?;
                pseudo.set_edgeuse_fillet_geom(primary, Some(solver_1.first_fillet_geom()?))?;
                self.edges.push(rib);

                // the 3-sided face's third edge, connecting the two rails in the
                // common face
                let across = pseudo.make_fillet_edge(Some(face_fv), Some(face_fv_2), Some(corner))?;
                pseudo
                    .fillet_edge_mut(across)?
                    .set_fillet_edge_type(FilletEdgeType::RailRailInterpolation);
                pseudo.fillet_edge_mut(across)?.set_original_face(common_face);
                self.edges.push(across);

                // a straight-across rib to end filletGeom2
                let rib = pseudo.make_fillet_edge(Some(face_fv_2), Some(side_fv), Some(corner))?;
                pseudo
                    .fillet_edge_mut(rib)?
                    .set_fillet_edge_type(FilletEdgeType::CrossSection);
                let primary = need(
                    eu1::primary_edgeuse(pseudo.brep(), rib)?,
                    "rib without primary edge use",
                )?;
                pseudo.set_edgeuse_fillet_geom(primary, Some(solver_2.first_fillet_geom()?))?;
                self.edges.push(rib);

                // all done
                return Ok(());
            } // end mixed-concavity case branch.

            if is_tangent_case || (mixed_convexity && angle_deg <= 5.0) {
                self.blending = true;
                // Two filleted edgeuses are tangent to each other: the new face vertex
                // is straight across from the rail-side edge intersection vertex.
                pseudo
                    .fillet_vertex_mut(face_fv)?
                    .set_fillet_vertex_type(FilletVertexType::Mate);
                pseudo.fillet_vertex_mut(face_fv)?.set_mate(0, Some(side_fv));
                pseudo.fillet_vertex_mut(side_fv)?.set_mate(0, Some(face_fv));
            } else {
                // Regular BEVEL corner: the new face vertex is the rails' intersection.
                pseudo
                    .fillet_vertex_mut(face_fv)?
                    .set_fillet_vertex_type(FilletVertexType::RailXRail);
                self.bevel = true;
            }

            // Setup rails of fillet
            self.nx2_rail_edgeuse(pseudo, target, corner, geometry_1, solver_1, rail, face_fv)?;
            self.nx2_rail_edgeuse(pseudo, target, corner, geometry_2, solver_2, other_rail, face_fv)?;
        } else if total == 4 {
            if is_tangent_case || (mixed_convexity && angle_deg <= 5.0) {
                // Two filleted edgeuses are tangent to each other
                self.blending = true;
                self.bevel = false;
            } else {
                // Regular BEVEL corner
                self.bevel = true;
                self.blending = false;
            }
        }

        // At this point, we have only one fillet (cross) edge at this corner.
        if self.vertices.len() != 2 {
            return Err(sm_err());
        }

        // Create the new cross-edge.
        let start = self.vertices.first().copied();
        let end = self.vertices.get(1).copied();
        let cross = pseudo.make_fillet_edge(start, end, Some(corner))?;
        let primary = need(
            eu1::primary_edgeuse(pseudo.brep(), cross)?,
            "cross edge without primary edge use",
        )?;
        let mate = need(
            eu1::mate(pseudo.brep(), primary)?,
            "cross edge use without mate",
        )?;
        pseudo.set_edgeuse_fillet_geom(primary, Some(geometry_1))?;
        pseudo.set_edgeuse_fillet_geom(mate, Some(geometry_2))?;

        // ('Or' is correct here: in the 3x2 case they are the same, not in 4x2.)
        let edge_type = if is_tangent_case || self.blending {
            FilletEdgeType::CrossSection
        } else {
            FilletEdgeType::FilletXFillet
        };
        pseudo.fillet_edge_mut(cross)?.set_fillet_edge_type(edge_type);
        self.edges.push(cross);

        // all done
        Ok(())
    }

    /// `SmFilletNx2Corner::Make4x2ConcaveCornerTopology`: the 3x2 mixed-convexity
    /// topology (a three-sided corner face) with one more, usually topological, edge at
    /// the corner vertex.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)] // one upstream routine
    fn make_4x2_concave_corner_topology(
        &mut self,
        k: &Kernel,
        target: &mut Brep,
        pseudo: &mut FilletBrep,
        corner: FilletCornerId,
        solvers: [&FilletSolver; 2],
        all: &[EdgeuseId],
        base_1: EdgeuseId,
        base_2: EdgeuseId,
        vertex: VertexId,
    ) -> Result<()> {
        // Get Solvers
        let [solver_1, solver_2] = solvers;
        let geometry_1 = solver_1.first_fillet_geom()?;
        let geometry_2 = solver_2.first_fillet_geom()?;

        let mut common_face: Option<FaceId> = None;
        let mut fv_on_common_face: Option<VertexId> = None;
        let mut fv_on_convex_edge: Option<VertexId> = None;
        let mut convex_geometry: Option<FilletGeomId> = None;
        let mut concave_geometry: Option<FilletGeomId> = None;

        // Make topology for all (both) of the non-filleted edges at this vertex.
        for ii in 0..4 {
            let current = *need(all.get(ii), "corner with fewer than four edge uses")?;
            if current == base_1 || current == base_2 {
                continue; // Do only the non-filleted ones.
            }
            // This edgeuse is in between two filleted edgeuses
            let mut fv = new_corner_vertex(k, pseudo, corner)?;
            self.vertices.push(fv);
            let side = need(
                find_side_edgeuse(target, base_1, current)?,
                "no side edge use beside the first filleted edge",
            )?;
            pseudo
                .fillet_vertex_mut(fv)?
                .set_point_class(PointClassType::Edgeuse, Some(EntityRef::Edgeuse(side)));
            pseudo
                .fillet_vertex_mut(fv)?
                .set_fillet_vertex_type(FilletVertexType::RailXEdgeuse);
            let mut fv_2 = fv;
            if check_face_complimentary_corner_concavity(
                k,
                target,
                vertex,
                base_1,
                base_2,
                Some(current),
            )? {
                let angle_1 = find_angle_between_edgeuses(k, target, vertex, current, base_1)?;
                let angle_2 = find_angle_between_edgeuses(k, target, vertex, current, base_2)?;
                if angle_1 > angle_2 {
                    common_face = face_by_two_edgeuses(target, Some(current), Some(base_1), None)?;
                    pseudo.fillet_vertex_mut(fv)?.set_point_class(
                        PointClassType::Edgeuse,
                        Some(EntityRef::Edgeuse(radial(target, side)?)),
                    );
                    let on_face = new_corner_vertex(k, pseudo, corner)?;
                    pseudo
                        .fillet_vertex_mut(on_face)?
                        .set_point_class(PointClassType::Face, common_face.map(EntityRef::Face));
                    fv_on_common_face = Some(on_face);
                    fv_2 = on_face;
                    pseudo
                        .fillet_vertex_mut(fv_2)?
                        .set_fillet_vertex_type(FilletVertexType::Mate);
                    core::mem::swap(&mut fv, &mut fv_2);
                    concave_geometry = Some(geometry_1);
                    convex_geometry = Some(geometry_2);
                } else {
                    common_face = face_by_two_edgeuses(target, Some(current), Some(base_2), None)?;
                    let on_face = new_corner_vertex(k, pseudo, corner)?;
                    pseudo
                        .fillet_vertex_mut(on_face)?
                        .set_point_class(PointClassType::Face, common_face.map(EntityRef::Face));
                    fv_on_common_face = Some(on_face);
                    fv_2 = on_face;
                    pseudo
                        .fillet_vertex_mut(fv_2)?
                        .set_fillet_vertex_type(FilletVertexType::Mate);
                    concave_geometry = Some(geometry_2);
                    convex_geometry = Some(geometry_1);
                }
            } else {
                fv_on_convex_edge = Some(fv);
            }
            let rail = solver_1.find_index_of_rail_x_side_edgeuse(target, side)?;
            self.nx2_rail_edgeuse(pseudo, target, corner, geometry_1, solver_1, rail, fv)?;
            let other_rail =
                solver_2.find_index_of_rail_x_side_edgeuse(target, radial(target, side)?)?;
            self.nx2_rail_edgeuse(pseudo, target, corner, geometry_2, solver_2, other_rail, fv_2)?;
        } // end loop on all (both) non-filleted edges.

        let on_face = fv_on_common_face.ok_or_else(|| SolidError::null_pointer(OP))?;
        let on_convex = fv_on_convex_edge.ok_or_else(|| SolidError::null_pointer(OP))?;
        pseudo.fillet_vertex_mut(on_face)?.set_mate(0, Some(on_convex));
        pseudo.fillet_vertex_mut(on_convex)?.set_mate(0, Some(on_face));
        self.vertices.push(on_face);

        // Will create a 3-sided corner face
        let mut prev = on_face; // i.e. m_vVertices[2]
        for ii in 0..3 {
            let current = *need(self.vertices.get(ii), "corner vertex index")?;
            let new_edge = pseudo.make_fillet_edge(Some(prev), Some(current), Some(corner))?;
            if ii == 1 || current == on_convex || prev == on_convex {
                pseudo
                    .fillet_edge_mut(new_edge)?
                    .set_fillet_edge_type(FilletEdgeType::CrossSection);
                let primary = need(
                    eu1::primary_edgeuse(pseudo.brep(), new_edge)?,
                    "fillet edge without primary edge use",
                )?;
                let geometry = if ii == 1 {
                    convex_geometry
                } else {
                    concave_geometry
                };
                pseudo.set_edgeuse_fillet_geom(primary, geometry)?;
            } else {
                pseudo
                    .fillet_edge_mut(new_edge)?
                    .set_fillet_edge_type(FilletEdgeType::RailRailInterpolation);
                pseudo.fillet_edge_mut(new_edge)?.set_original_face(common_face);
            }
            self.edges.push(new_edge);
            prev = current;
        }
        Ok(())
    }
}
