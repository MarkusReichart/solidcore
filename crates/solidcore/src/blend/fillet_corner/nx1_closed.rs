//! Topology at the vertex of a closed filleted edge (the Nx1 closed corner).
//!
//! Upstream: `source/SMLib/src/SmFilletCorner.cpp` (4962–5129).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`FilletCorner::make_nx1_closed_corner_topology`] | `SmFilletNx1ClosedCorner::MakeCornerTopology` |
//! | [`Nx1ClosedTopologyEnv`] | the executive's pseudo-B-rep and `GetFilletSolverOfEdgeuse` |
//!
//! Deviations: the caller lends the target, the pseudo-B-rep and the loaded solver of
//! the filleted edge use. Missing topology that upstream dereferences is an input
//! error. There is no recursion.

use crate::blend::fillet_corner::p1::test_g1_continuity;
use crate::blend::fillet_corner::p3::{FilletCorner, FilletCornerKind};
use crate::blend::fillet_geom::p1::{
    CornerLookup, FilletBrep, FilletCornerId, FilletEdgeType, FilletVertex, FilletVertexType,
};
use crate::blend::fillet_solver::p1::FilletSolver;
use crate::error::{EntityRef, Result, SolidError};
use crate::intersect::point_class::PointClassType;
use crate::kernel::Kernel;
use crate::math::{rad_to_deg, Point3};
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::topology::edge;
use crate::topology::edgeuse::{p1 as eu1, p3 as eu3};
use crate::topology::{Brep, EdgeId, EdgeuseId, Vertex, VertexId};

const OP: &str = "blend::fillet_corner::nx1_closed::make_nx1_closed_corner_topology";

/// What `SmFilletNx1ClosedCorner::MakeCornerTopology` works on.
#[derive(Debug)]
pub struct Nx1ClosedTopologyEnv<'a> {
    /// The executive's kernel.
    pub k: &'a Kernel,
    /// The B-rep being filleted (`EvaluateBinormal` may build a face's surface cache
    /// on it).
    pub target: &'a mut Brep,
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

fn edge_of(target: &Brep, eu: EdgeuseId) -> Result<EdgeId> {
    need(eu1::edge(target, eu)?, "edge use without edge")
}

/// The parameter at the start of the edge use: its edge interval's minimum when it runs
/// with the edge, else the maximum.
fn start_parameter(target: &Brep, eu: EdgeuseId) -> Result<f64> {
    let ivl = edge::interval(target, edge_of(target, eu)?)?;
    Ok(if eu1::orientation(target, eu)? == OrientType::Same {
        ivl.min()
    } else {
        ivl.max()
    })
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

/// The corner store lookup that resolves this corner (`this`) to `c`.
fn lookup(corner: FilletCornerId, c: &FilletCorner) -> CornerLookup<'_> {
    CornerLookup {
        this: Some((Some(corner), c)),
        ..Default::default()
    }
}

impl FilletCorner {
    /// `SmFilletNx1ClosedCorner::MakeCornerTopology`: the corner at the vertex of a
    /// closed filleted edge. Every other edge use at the vertex (a seam), taken on the
    /// side the filleted edge's start tangent points into, gets a rail-by-edge-use
    /// vertex where both senses of the rail crossing it start. With two edge uses at
    /// the vertex a second vertex goes on the filleted edge use's face (or its radial's
    /// when that is the seam's): the first one's mate when the edge is G1 closed, else
    /// the fillet/fillet intersection, both senses of its rail starting there. The two
    /// first vertices are joined by a cross-section edge (G1, blending) or a
    /// fillet/fillet edge (the fillet geometry is then split), both its uses on the
    /// solver's first fillet geometry.
    ///
    /// # Errors
    ///
    /// Invalid topology or a different corner kind; `SM_ERR` when the fillet edge does
    /// not have two edge uses; those of `sm_TestG1Continuity`, `EvaluateBinormal`, the
    /// curve evaluation and the pseudo-topology operations.
    pub fn make_nx1_closed_corner_topology(&mut self, env: Nx1ClosedTopologyEnv<'_>) -> Result<()> {
        let Nx1ClosedTopologyEnv {
            k,
            target,
            pseudo,
            corner,
            solver,
        } = env;
        if self.kind != FilletCornerKind::NByOneClosed {
            return Err(invalid("requires an Nx1 closed corner"));
        }
        let filleted = need(
            self.solver_eus
                .as_deref()
                .and_then(|eus| eus.first())
                .copied(),
            "corner without filleted edge use",
        )?;
        let geometry = solver.first_fillet_geom()?;
        let all = need(self.all_eus.clone(), "corner without outward edge uses")?;
        let total = all.len();
        let filleted_edge = edge_of(target, filleted)?;

        self.g1_case = test_g1_continuity(
            k,
            target,
            self.vertex,
            Some(filleted_edge),
            Some(filleted_edge),
            rad_to_deg(self.tangency_tol_radians),
        )?
        .0;


        // for every edge at this vertex, except the filleted edge
        let mut seam = None;
        for eu in all {
            if eu == filleted {
                continue;
            }
            // the seam edge use on the side the filleted edge's start tangent points into
            let t = start_parameter(target, eu)?;
            let binormal = eu3::evaluate_binormal(
                k,
                target,
                eu,
                t,
                false,
                &eu3::BinormalOptions::default(),
            )?
            .binormal;
            let t1 = start_parameter(target, filleted)?;
            let curve = need(edge::curve(target, filleted_edge)?, "edge without curve")?;
            let pv = curve.evaluate(k, t1, 1, true, true)?;
            let tangent = need(pv.get(1).copied(), "curve without first derivative")?;
            let mut seam_eu = eu;
            if tangent.dot(binormal) < 0.0 {
                seam_eu = need(eu1::radial(target, eu)?, "seam edge use without radial")?;
            }
            seam = Some(seam_eu);

            let fv = new_corner_vertex(k, pseudo, corner)?;
            let record = pseudo.fillet_vertex_mut(fv)?;
            record.set_point_class(PointClassType::Edgeuse, Some(EntityRef::Edgeuse(seam_eu)));
            record.set_fillet_vertex_type(FilletVertexType::RailXEdgeuse);
            self.vertices.push(fv);

            let rail = solver.find_index_of_rail_x_side_edgeuse(target, seam_eu)?;
            for orient in [OrientType::Same, OrientType::Opposite] {
                pseudo.make_rail_edgeuse_on(
                    geometry,
                    rail,
                    Some(fv),
                    orient,
                    None,
                    lookup(corner, self),
                    solver,
                    target,
                )?;
            }
        }

        let fv0 = need(self.vertices.first().copied(), "corner without a seam vertex")?;
        if total == 2 {
            let fv = new_corner_vertex(k, pseudo, corner)?;
            let seam = need(seam, "corner without a seam edge use")?;
            let mut face = eu1::face(target, filleted)?;
            if face == eu1::face(target, seam)? {
                let radial = need(eu1::radial(target, filleted)?, "edge use without radial")?;
                face = eu1::face(target, radial)?;
            }
            pseudo
                .fillet_vertex_mut(fv)?
                .set_point_class(PointClassType::Face, face.map(EntityRef::Face));
            if self.g1_case {
                let record = pseudo.fillet_vertex_mut(fv)?;
                record.set_fillet_vertex_type(FilletVertexType::Mate);
                record.set_mate(0, Some(fv0));
                pseudo.fillet_vertex_mut(fv0)?.set_mate(0, Some(fv));
            } else {
                // for non-g1 closed curves the vertex is at the fillet/fillet
                // intersection
                pseudo
                    .fillet_vertex_mut(fv)?
                    .set_fillet_vertex_type(FilletVertexType::FilletXFillet);
            }
            self.vertices.push(fv);

            let rail = solver.find_index_of_rail_on_face(target, face)?;
            for orient in [OrientType::Same, OrientType::Opposite] {
                pseudo.make_rail_edgeuse_on(
                    geometry,
                    rail,
                    Some(fv),
                    orient,
                    None,
                    lookup(corner, self),
                    solver,
                    target,
                )?;
            }
        }

        let fv1 = self.vertices.get(1).copied();
        let edge = pseudo.make_fillet_edge(Some(fv0), fv1, Some(corner))?;
        if self.g1_case {
            pseudo
                .fillet_edge_mut(edge)?
                .set_fillet_edge_type(FilletEdgeType::CrossSection);
            self.blending = true;
        } else {
            // a closed curve not tangent: the edge's shape is the fillet/fillet
            // intersection
            pseudo.geom_mut(geometry)?.need_split = true;
            pseudo
                .fillet_edge_mut(edge)?
                .set_fillet_edge_type(FilletEdgeType::FilletXFillet);
            self.blending = false;
        }

        if pseudo.brep().edge(edge)?.edgeuses.len() != 2 {
            return Err(SolidError::sm_err(OP));
        }
        let primary = need(
            eu1::primary_edgeuse(pseudo.brep(), edge)?,
            "fillet edge without primary edge use",
        )?;
        pseudo.set_edgeuse_fillet_geom(primary, Some(geometry))?;
        let mate = need(eu1::mate(pseudo.brep(), primary)?, "fillet edge use without mate")?;
        pseudo.set_edgeuse_fillet_geom(mate, Some(geometry))?;
        self.edges.push(edge);
        Ok(())
    }
}
