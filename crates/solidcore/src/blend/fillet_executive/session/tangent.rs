//! The live `ProcessSmoothEdgeSequences` access of a fillet session: the queries of
//! `CollectTangentEdgeSequence` and the edits of `CreateSmoothEdgeSequence` over the
//! session's own solvers, corners and target B-rep.
//!
//! Upstream: `source/SMLib/src/SmFilletExecutive.cpp` (`ProcessSmoothEdgeSequences`,
//! `CollectTangentEdgeSequence`, `CreateSmoothEdgeSequence`),
//! `source/SMLib/src/SmFilletSolver.cpp` (`GetFilletCorner`, `GetOtherFilletCorner`),
//! `source/SMLib/src/SmFilletCorner.cpp` (`GetTangentFilletSolver`).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | `SessionTangent` (crate) | the executive's solvers, corners and filleted edges as `ProcessSmoothEdgeSequences` reaches them |
//! | `SessionTangent::new` (crate) | (the session's executive maps in the corner code's form) |
//!
//! Deviations: the corner lookups read a copy of the executive's vertex map taken when
//! the access is made (`ProcessSmoothEdgeSequences` changes no map), and
//! `GetTangentFilletSolver` runs on the corner code's [`FilletExecutiveData`], of which
//! only the edge-use solver map is filled (the only member it reads). The edge's
//! `SetCurve(.., TRUE)` notifies the session's target records without a context (no
//! listener sees it) and drops the session's curve cache of that edge, which the new
//! curve replaces. Upstream swaps a curve-based solver's original curve only when it is
//! the edge's curve object (a pointer comparison); the session's solvers hold their own
//! copy, set from the edge's curve, so a set original curve is swapped. No recursion.

use std::collections::BTreeMap;

use super::{invalid, FilletBody, FilletObject, SessionObjects, SolverObject};
use crate::blend::fillet_corner::p3::FilletExecutiveData;
use crate::blend::fillet_executive::p2::{SmoothSequenceAccess, TangentSequenceAccess};
use crate::blend::fillet_executive::p3::{FilletExecutive, VertexEntry};
use crate::blend::fillet_solver::p1::FilletSolver;
use crate::curve::curve::Curve;
use crate::error::{EntityRef, Result, SolidError};
use crate::geom::extent1d::Extent1d;
use crate::kernel::Kernel;
use crate::topology::{edge, edgeuse::p1 as eu1, Brep, BrepId, EdgeId};

const OP: &str = "blend::fillet_executive::session::tangent";

/// The session's solvers, corners and target as `ProcessSmoothEdgeSequences` reads and
/// edits them.
pub(super) struct SessionTangent<'a, 'target> {
    objects: &'a mut SessionObjects<'target>,
    k: &'a Kernel,
    vertex_corners: BTreeMap<Option<FilletObject>, VertexEntry<FilletObject>>,
    data: FilletExecutiveData,
}

impl<'a, 'target> SessionTangent<'a, 'target> {
    /// The access over `objects`, with `executive`'s vertex map and edge-use solver map.
    ///
    /// # Errors
    /// An edge-use solver map entry that is not a target edge use and a session solver.
    pub(super) fn new(
        objects: &'a mut SessionObjects<'target>,
        k: &'a Kernel,
        executive: &FilletExecutive<FilletObject>,
    ) -> Result<Self> {
        let mut edgeuse_solver_map = BTreeMap::new();
        for (&eu, &solver) in &executive.edgeuse_solver_map {
            let FilletObject::Topology(FilletBody::Target, EntityRef::Edgeuse(eu)) = eu else {
                return Err(invalid("edge-use solver map entry is not a target edge use"));
            };
            let FilletObject::Solver(solver) = solver else {
                return Err(invalid("edge-use solver map entry is not a session solver"));
            };
            edgeuse_solver_map.insert(eu, solver);
        }
        Ok(Self {
            objects,
            k,
            vertex_corners: executive.vertex_corner_map.clone(),
            data: FilletExecutiveData {
                edgeuse_solver_map,
                ..FilletExecutiveData::default()
            },
        })
    }

    fn solver(&self, solver: FilletObject) -> Result<&FilletSolver> {
        let FilletObject::Solver(id) = solver else {
            return Err(invalid("tangent object is not an actual solver"));
        };
        Ok(self.objects.solver(id)?.base())
    }

    /// `pFS->GetEdgeuse(0)->GetEdge()`.
    fn solver_edge(&self, solver: FilletObject) -> Result<EdgeId> {
        let eu = self
            .solver(solver)?
            .edgeuse(0)?
            .ok_or_else(|| SolidError::null_pointer(OP))?;
        eu1::edge(self.objects.target, eu)?.ok_or_else(|| SolidError::null_pointer(OP))
    }
}

impl TangentSequenceAccess<FilletObject> for SessionTangent<'_, '_> {
    fn tangent_solver_edge(&self, solver: FilletObject) -> Result<FilletObject> {
        let e = self.solver_edge(solver)?;
        Ok(FilletObject::Topology(FilletBody::Target, EntityRef::Edge(e)))
    }

    fn tangent_edge_edgeuses(&self, edge_object: FilletObject) -> Result<Vec<FilletObject>> {
        let FilletObject::Topology(FilletBody::Target, EntityRef::Edge(e)) = edge_object else {
            return Err(invalid("tangent object is not a target edge"));
        };
        Ok(edge::edgeuses(self.objects.target, e, None, None)?
            .into_iter()
            .map(|eu| FilletObject::Topology(FilletBody::Target, EntityRef::Edgeuse(eu)))
            .collect())
    }

    fn tangent_solver_corner(
        &self,
        solver: FilletObject,
        index: usize,
    ) -> Result<Option<FilletObject>> {
        if index > 1 {
            return Ok(None);
        }
        let vertex = self.solver(solver)?.vertex(self.objects.target, index)?;
        let key = vertex.map(|v| FilletObject::Topology(FilletBody::Target, EntityRef::Vertex(v)));
        // GetFilletCornerOfVertex
        match self.vertex_corners.get(&key) {
            None => Ok(None),
            Some(VertexEntry::Corner(corner)) => Ok(Some(*corner)),
            Some(VertexEntry::SolverEdgeuses(_)) => Err(SolidError::invalid_input(
                OP,
                "vertex still contains a pending solver list",
            )),
        }
    }

    fn tangent_solver_other_corner(
        &self,
        solver: FilletObject,
        corner: FilletObject,
    ) -> Result<Option<FilletObject>> {
        let corner0 = self.tangent_solver_corner(solver, 0)?;
        let corner1 = self.tangent_solver_corner(solver, 1)?;
        if Some(corner) == corner0 {
            return Ok(corner1);
        }
        if Some(corner) == corner1 {
            return Ok(corner0);
        }
        // Error: pThisCorner is not one of ours.
        Ok(None)
    }

    fn tangent_corner_tangent_solver(
        &self,
        corner: FilletObject,
        solver: FilletObject,
    ) -> Result<Option<FilletObject>> {
        let (FilletObject::Corner(cid), FilletObject::Solver(sid)) = (corner, solver) else {
            return Err(invalid("tangent objects are not an actual corner and solver"));
        };
        let corner = self
            .objects
            .corners
            .get(cid.0 as usize)
            .ok_or_else(|| invalid("corner is not live in this session"))?;
        Ok(corner
            .tangent_fillet_solver(self.k, self.objects.target, &self.data, sid)?
            .map(FilletObject::Solver))
    }
}

impl SmoothSequenceAccess<FilletObject> for SessionTangent<'_, '_> {
    fn smooth_edge(&self, solver: FilletObject) -> Result<EdgeId> {
        self.solver_edge(solver)
    }

    fn smooth_target(&mut self) -> &mut Brep {
        &mut *self.objects.target
    }

    fn smooth_set_edge_curve(&mut self, _k: &Kernel, e: EdgeId, curve: Curve) -> Result<()> {
        let serial = self.objects.target_serial;
        let key_of = |entity| crate::topology::brep::p3::topology_ref(serial, entity);
        crate::topology::brep::p3::seed_records(
            &mut self.objects.notification_records,
            self.objects.target,
            &key_of,
        );
        let mut env = edge::NotifyEnv {
            records: &mut self.objects.notification_records,
            attributes: &mut self.objects.attributes,
            context: None,
            key_of: &key_of,
        };
        edge::set_curve_with_body_identity(
            self.objects.target,
            e,
            Some(curve),
            edge::SetCurveOptions {
                delete_old_curve: true,
                ..edge::SetCurveOptions::default()
            },
            &mut env,
            [None, None],
            key_of(EntityRef::Brep(BrepId::from_index(0))),
        )?;
        self.objects.edge_caches.remove(&e);
        Ok(())
    }

    fn smooth_replace_original_curve(
        &mut self,
        solver: FilletObject,
        e: EdgeId,
        curve: &Curve,
    ) -> Result<()> {
        let FilletObject::Solver(id) = solver else {
            return Err(invalid("tangent object is not an actual solver"));
        };
        let _ = e;
        // SM_CAST_PTR(SmCurveBasedFS, pFS): the variable-radius solver is the curve-based
        // one the session keeps; its original curve, when set, is its edge's (see the
        // module's Deviations)
        if let SolverObject::VariableRadius(fs) = self.objects.solver_mut(id)? {
            if fs.curve_based.original_curve().is_some() {
                fs.curve_based.set_original_curve(Some(curve.clone()));
            }
        }
        Ok(())
    }

    fn smooth_set_law_interval(&mut self, solver: FilletObject, interval: Extent1d) -> Result<()> {
        let FilletObject::Solver(id) = solver else {
            return Err(invalid("tangent object is not an actual solver"));
        };
        if let SolverObject::VariableRadius(fs) = self.objects.solver_mut(id)? {
            fs.set_fillet_law_interval(interval);
        }
        Ok(())
    }
}
