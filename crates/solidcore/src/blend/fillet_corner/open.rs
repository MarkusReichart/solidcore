//! Topology of the open fillet corner (an end at a lamina edge).
//!
//! Upstream: `source/SMLib/src/SmFilletCorner.cpp` (8699–8748).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`FilletCorner::make_open_corner_topology`] | `SmFilletOpenCorner::MakeCornerTopology` |
//! | [`OpenTopologyEnv`] | the executive's pseudo-B-rep and `GetFilletSolverOfEdgeuse` |
//!
//! Deviations: the caller lends the target, the pseudo-B-rep and the loaded solver of
//! the filleted edge use. Missing topology that upstream dereferences is an input
//! error. There is no recursion.

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
use crate::topology::{Brep, Vertex};

const OP: &str = "blend::fillet_corner::open::make_open_corner_topology";

/// What `SmFilletOpenCorner::MakeCornerTopology` works on.
#[derive(Debug)]
pub struct OpenTopologyEnv<'a> {
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
    /// `SmFilletOpenCorner::MakeCornerTopology`: the end of a fillet at a lamina edge
    /// (shell filleting). The filleted edge use's clockwise neighbour and its radial's
    /// counter-clockwise neighbour each get a vertex where the fillet's rail crosses
    /// them, joined by the fillet-end edge on the solver's first fillet geometry.
    ///
    /// # Errors
    ///
    /// Invalid topology or a different corner kind; those of the pseudo-topology
    /// operations.
    pub fn make_open_corner_topology(&mut self, env: OpenTopologyEnv<'_>) -> Result<()> {
        let OpenTopologyEnv {
            k,
            target,
            pseudo,
            corner,
            solver,
        } = env;
        if self.kind != FilletCornerKind::Open {
            return Err(invalid("requires an open corner"));
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
        let geometry = solver.first_fillet_geom()?;
        for eu in adjacent {
            let fv = pseudo.brep_mut().add_unlinked_vertex_record(Vertex {
                point: Point3::ZERO,
                tolerance: 0.0,
                ..Vertex::default()
            })?;
            pseudo.insert_fillet_vertex(fv, FilletVertex::new(k, Some(corner)));
            let record = pseudo.fillet_vertex_mut(fv)?;
            record.set_point_class(PointClassType::Edgeuse, Some(EntityRef::Edgeuse(eu)));
            record.set_fillet_vertex_type(FilletVertexType::RailXEdgeuse);
            self.vertices.push(fv);
            let rail = solver.find_index_of_rail_x_side_edgeuse(target, eu)?;
            // Create rail-edgeuse using this vertex
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
        }

        let start = self.vertices.first().copied();
        let end = self.vertices.get(1).copied();
        let edge = pseudo.make_fillet_edge(start, end, Some(corner))?;
        pseudo
            .fillet_edge_mut(edge)?
            .set_fillet_edge_type(FilletEdgeType::FilletEnd);
        let primary = need(
            eu1::primary_edgeuse(pseudo.brep(), edge)?,
            "fillet edge without primary edge use",
        )?;
        pseudo.set_edgeuse_fillet_geom(primary, Some(geometry))?;
        self.edges.push(edge);
        Ok(())
    }
}
