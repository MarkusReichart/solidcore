//! Topology of the tangent 1x1 and 2x2 fillet corners: one cross-section edge.
//!
//! Upstream: `source/SMLib/src/SmFilletCorner.cpp` (4372–4436, 4494–4561).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`FilletCorner::make_1x1_corner_topology`] | `SmFillet1x1Corner::MakeCornerTopology` |
//! | [`FilletCorner::make_2x2_corner_topology`] | `SmFillet2x2Corner::MakeCornerTopology` |
//! | [`TangentTopologyEnv`] | the executive's pseudo-B-rep and `GetFilletSolverOfEdgeuse` |
//!
//! Deviations: the caller lends the target, the pseudo-B-rep and the loaded solvers of
//! the corner's filleted edge uses (upstream asks the executive; a missing solver is
//! `NER` there and an input error here). Missing topology that upstream dereferences is
//! an input error. There is no recursion.

use crate::blend::fillet_corner::p3::{FilletCorner, FilletCornerKind};
use crate::blend::fillet_geom::p1::{
    CornerLookup, FilletBrep, FilletCornerId, FilletEdgeType, FilletVertex, FilletVertexType,
};
use crate::blend::fillet_solver::p1::FilletSolver;
use crate::error::{EntityRef, Result, SolidError};
use crate::intersect::point_class::PointClassType;
use crate::kernel::Kernel;
use crate::math::Point3;
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::topology::edgeuse::p1 as eu1;
use crate::topology::{Brep, EdgeuseId, FaceId, Vertex, VertexId};

const OP: &str = "blend::fillet_corner::one_by_one";

/// What `SmFillet1x1Corner::MakeCornerTopology` and `SmFillet2x2Corner::MakeCornerTopology`
/// work on.
#[derive(Debug)]
pub struct TangentTopologyEnv<'a> {
    /// The executive's kernel.
    pub k: &'a Kernel,
    /// The B-rep being filleted (the corner's edge uses and faces).
    pub target: &'a Brep,
    /// The executive's pseudo-B-rep.
    pub pseudo: &'a mut FilletBrep,
    /// This corner's identity in the executive's corner store.
    pub corner: FilletCornerId,
    /// The solver of the first filleted edge use (`GetFilletSolverOfEdgeuse`).
    pub solver: &'a FilletSolver,
    /// The solver of the second filleted edge use (a 2x2 corner's).
    pub second_solver: Option<&'a FilletSolver>,
}

fn invalid(reason: &'static str) -> SolidError {
    SolidError::invalid_input(OP, reason)
}

fn need<T>(value: Option<T>, reason: &'static str) -> Result<T> {
    value.ok_or_else(|| invalid(reason))
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

/// The filleted edge use `index` of the corner (`(*m_pSolverEUs)[index]`).
fn solver_eu(corner: &FilletCorner, index: usize) -> Result<EdgeuseId> {
    need(
        corner
            .solver_eus
            .as_deref()
            .and_then(|eus| eus.get(index))
            .copied(),
        "corner without its filleted edge use",
    )
}

impl FilletCorner {
    /// The two fillet vertices on either side of the filleted edge use `filleted` (on
    /// its face, a cross-section vertex, and on its radial's face, its mate) and the
    /// cross-section fillet edge joining them: the part `SmFillet1x1Corner` and
    /// `SmFillet2x2Corner::MakeCornerTopology` share. Returns the two vertices, the
    /// first one's face and the edge.
    fn make_cross_section_corner(
        &mut self,
        env: &mut TangentTopologyEnv<'_>,
        filleted: EdgeuseId,
    ) -> Result<([VertexId; 2], Option<FaceId>, crate::topology::EdgeId)> {
        let (k, target, corner) = (env.k, env.target, env.corner);

        // Create a fillet vertex on one side of the edge
        let fv0 = new_corner_vertex(k, env.pseudo, corner)?;
        let face0 = eu1::face(target, filleted)?;
        let fv = env.pseudo.fillet_vertex_mut(fv0)?;
        fv.set_point_class(PointClassType::Face, face0.map(EntityRef::Face));
        fv.set_fillet_vertex_type(FilletVertexType::OnCrossSection);
        self.vertices.push(fv0);

        // Create another fillet vertex on the other side of the edge
        let fv1 = new_corner_vertex(k, env.pseudo, corner)?;
        let radial = need(eu1::radial(target, filleted)?, "edge use without radial")?;
        let face1 = eu1::face(target, radial)?;
        let fv = env.pseudo.fillet_vertex_mut(fv1)?;
        fv.set_point_class(PointClassType::Face, face1.map(EntityRef::Face));
        fv.set_fillet_vertex_type(FilletVertexType::Mate);
        fv.set_mate(0, Some(fv0));
        env.pseudo.fillet_vertex_mut(fv0)?.set_mate(0, Some(fv1));
        self.vertices.push(fv1);

        // Corner consists of one fillet edge only (joining pFV0 & pFV1)
        let edge = env.pseudo.make_fillet_edge(Some(fv0), Some(fv1), Some(corner))?;
        env.pseudo
            .fillet_edge_mut(edge)?
            .set_fillet_edge_type(FilletEdgeType::CrossSection);
        Ok(([fv0, fv1], face0, edge))
    }

    /// The primary edge use of `edge` and its mate (pseudo-B-rep).
    fn edge_uses(pseudo: &FilletBrep, edge: crate::topology::EdgeId) -> Result<[EdgeuseId; 2]> {
        let primary = need(
            eu1::primary_edgeuse(pseudo.brep(), edge)?,
            "fillet edge without primary edge use",
        )?;
        let mate = need(eu1::mate(pseudo.brep(), primary)?, "fillet edge use without mate")?;
        Ok([primary, mate])
    }

    /// `SmFillet1x1Corner::MakeCornerTopology`: a closed edge separating two adjacent
    /// faces at a vertex of its own (two tori meeting in an '8'). One fillet edge joins
    /// a cross-section vertex on the filleted edge use's face to its mate on the
    /// radial's face, each fillet edge use on the solver's first fillet geometry, and
    /// each vertex starts both senses of its (closed) rail.
    ///
    /// # Errors
    ///
    /// Invalid topology or a different corner kind; `SM_ERR` when the fillet edge does
    /// not have two edge uses; those of the pseudo-topology operations.
    pub fn make_1x1_corner_topology(&mut self, mut env: TangentTopologyEnv<'_>) -> Result<()> {
        if self.kind != FilletCornerKind::OneByOne {
            return Err(invalid("requires a 1x1 corner"));
        }
        let filleted = solver_eu(self, 0)?;
        let ([fv0, fv1], face0, edge) = self.make_cross_section_corner(&mut env, filleted)?;
        let TangentTopologyEnv {
            target,
            pseudo,
            corner,
            solver,
            ..
        } = env;

        if pseudo.brep().edge(edge)?.edgeuses.len() != 2 {
            return Err(SolidError::sm_err(OP));
        }
        let geometry = solver.first_fillet_geom()?;
        for eu in Self::edge_uses(pseudo, edge)? {
            pseudo.set_edgeuse_fillet_geom(eu, Some(geometry))?;
        }
        self.edges.push(edge);

        // Setup rails of fillet: two edge uses for each (closed) rail
        let rail = solver.find_index_of_rail_on_face(target, face0)?;
        for (index, fv) in [(rail, fv0), (1 - rail, fv1)] {
            for orient in [OrientType::Same, OrientType::Opposite] {
                pseudo.make_rail_edgeuse_on(
                    geometry,
                    index,
                    Some(fv),
                    orient,
                    None,
                    CornerLookup {
                        this: Some((Some(corner), self)),
                        ..Default::default()
                    },
                    solver,
                    target,
                )?;
            }
        }

        self.blending = true;
        Ok(())
    }

    /// `SmFillet2x2Corner::MakeCornerTopology`: two filleted edges meeting with G1
    /// continuity at the corner. One cross-section fillet edge separates the two
    /// fillets: its primary edge use on the first solver's first fillet geometry, its
    /// mate on the second's, and each fillet geometry's rails end at the two vertices.
    ///
    /// # Errors
    ///
    /// Invalid topology or a different corner kind; a missing second solver (`NER`
    /// upstream); those of the pseudo-topology operations.
    pub fn make_2x2_corner_topology(&mut self, mut env: TangentTopologyEnv<'_>) -> Result<()> {
        if self.kind != FilletCornerKind::TwoByTwo {
            return Err(invalid("requires a 2x2 corner"));
        }
        let filleted = solver_eu(self, 0)?;
        solver_eu(self, 1)?;
        let second = need(env.second_solver, "2x2 corner without its second solver")?;
        let ([fv0, fv1], face0, edge) = self.make_cross_section_corner(&mut env, filleted)?;
        let TangentTopologyEnv {
            target,
            pseudo,
            corner,
            solver,
            ..
        } = env;

        let [primary, mate] = Self::edge_uses(pseudo, edge)?;
        let geometry1 = solver.first_fillet_geom()?;
        pseudo.set_edgeuse_fillet_geom(primary, Some(geometry1))?;
        let geometry2 = second.first_fillet_geom()?;
        pseudo.set_edgeuse_fillet_geom(mate, Some(geometry2))?;
        self.edges.push(edge);

        // Setup rails of fillet
        for (geometry, rail_solver) in [(geometry1, solver), (geometry2, second)] {
            let rail = rail_solver.find_index_of_rail_on_face(target, face0)?;
            for (index, fv) in [(rail, fv0), (1 - rail, fv1)] {
                pseudo.make_rail_edgeuse_on(
                    geometry,
                    index,
                    Some(fv),
                    OrientType::Unknown,
                    None,
                    CornerLookup {
                        this: Some((Some(corner), self)),
                        ..Default::default()
                    },
                    rail_solver,
                    target,
                )?;
            }
        }

        self.blending = true;
        Ok(())
    }
}
