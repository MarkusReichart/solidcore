//! The corner code's [`CornerEnv`] over a session's own objects.
//!
//! Upstream: `source/SMLib/src/SmFilletCorner.cpp`, `source/SMLib/src/SmFilletExecutive.cpp`
//! (the corners reach `m_pExecutive`'s members and solvers through pointers).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | `SessionObjects::with_corner_env` (crate) | a corner method running on the executive's live members |
//! | `SessionObjects::tsect_surfaces` (private) | `SmTsectPnt::m_apUserPointer[0..2]` of each rail vertex use (`SmFilletGeom.cpp:10079-10091`) |
//!
//! Deviations: the corner code ([`CornerEnv`], `blend.fillet_corner_*`) reads the
//! executive's members as [`FilletExecutiveData`] over typed ids and the solvers as
//! one slice of base [`FilletSolver`]s indexed by [`FilletSolverId`]; the session keeps
//! them as [`FilletExecutive<FilletObject>`] and concrete [`SolverObject`]s. For one
//! corner call the session's state moves into that form (the solvers' bases are
//! detached from their concrete solvers, the corner table, extended-surface map,
//! topology edges, global-merge switch, maps and solver list are converted) and moves
//! back afterwards, on success and failure alike: [`SessionObjects::with_corner_env`]
//! is the only way in and always restores. There is one corner store; the running
//! corner is out of the table while it runs (`FilletExecutiveData::corners` holds the
//! others). The intersection points' user surfaces are not modelled on
//! [`TsectPnt`](crate::intersect::surface_intersector::p1::TsectPnt): they are
//! recomputed per rail vertex use from the solver's edge uses, as
//! `MakeRailEdgeuseOn` set them. A corner error is copied to the executive when the
//! corner call recorded one. No recursion.

use std::collections::BTreeMap;

use super::{invalid, FilletBody, FilletObject, SessionObjects, SolverObject};
use crate::blend::fillet_corner::p3::{
    CornerEnv, ExecutiveErrorObject, FilletCorner, FilletExecutiveData,
};
use crate::blend::fillet_executive::p3::{ErrorInfo, FilletExecutive, VertexEntry};
use crate::blend::fillet_geom::p1::{FilletCornerId, FilletEdgeType};
use crate::blend::fillet_solver::p1::FilletSolver;
use crate::error::{EntityRef, Result};
use crate::kernel::Kernel;
use crate::surface::surface::Surface;
use crate::topology::{EdgeId, EdgeuseId, VertexId, VertexuseId, VertexuseOwner};

/// The session state moved into corner form for one corner call.
struct Taken {
    corner: FilletCorner,
    data: FilletExecutiveData,
    solvers: Vec<FilletSolver>,
}

fn target_edgeuse(object: FilletObject) -> Option<EdgeuseId> {
    match object {
        FilletObject::Topology(FilletBody::Target, EntityRef::Edgeuse(eu)) => Some(eu),
        _ => None,
    }
}

fn target_vertex(object: Option<FilletObject>) -> Option<VertexId> {
    match object {
        Some(FilletObject::Topology(FilletBody::Target, EntityRef::Vertex(v))) => Some(v),
        _ => None,
    }
}

fn pseudo_edge(object: FilletObject) -> Option<EdgeId> {
    match object {
        FilletObject::Topology(FilletBody::Pseudo, EntityRef::Edge(e)) => Some(e),
        _ => None,
    }
}

/// An error object of the corner code as the session names it (`None`: an object the
/// session has no handle for, the temporary end plane).
fn error_object(object: ExecutiveErrorObject) -> Option<FilletObject> {
    use super::insertion::geometry::InsertionBody;
    match object {
        ExecutiveErrorObject::TargetVertex(v) => Some(FilletObject::Topology(
            FilletBody::Target,
            EntityRef::Vertex(v),
        )),
        ExecutiveErrorObject::FilletVertex(v) => Some(FilletObject::Topology(
            FilletBody::Pseudo,
            EntityRef::Vertex(v),
        )),
        ExecutiveErrorObject::FilletGeomSurface(g) => Some(FilletObject::Geometry(g)),
        ExecutiveErrorObject::TargetFaceSurface { face, .. } => {
            Some(FilletObject::FaceSurface(InsertionBody::Target, face))
        }
        ExecutiveErrorObject::TempEndPlane(_) => None,
    }
}

fn error_info(source: &ErrorInfo<ExecutiveErrorObject>) -> ErrorInfo<FilletObject> {
    let one = |o: Option<ExecutiveErrorObject>| o.and_then(error_object);
    let many = |v: &[Option<ExecutiveErrorObject>]| -> Vec<Option<FilletObject>> {
        v.iter().map(|o| one(*o)).collect()
    };
    let [fv0, fv1] = source.fillet_vertices;
    ErrorInfo {
        code: source.code,
        status: source.status,
        corner_type: source.corner_type,
        edge_type: source.edge_type,
        executive: Some(FilletObject::Executive),
        vertex: one(source.vertex),
        edge: one(source.edge),
        side_face: one(source.side_face),
        end_face: one(source.end_face),
        fillet_vertices: [one(fv0), one(fv1)],
        surfaces: many(&source.surfaces),
        edges: many(&source.edges),
        curves: many(&source.curves),
        location: source.location,
        comment: source.comment.clone(),
        entity_name: source.entity_name.clone(),
    }
}

fn detach_bases(slots: &mut [Option<SolverObject>]) -> Vec<FilletSolver> {
    slots
        .iter_mut()
        .map(|slot| match slot {
            Some(solver) => {
                std::mem::replace(solver.base_mut(), FilletSolver::new(0.0, 0.0, 0.0))
            }
            None => FilletSolver::new(0.0, 0.0, 0.0),
        })
        .collect()
}

fn attach_bases(slots: &mut [Option<SolverObject>], bases: Vec<FilletSolver>) {
    for (slot, base) in slots.iter_mut().zip(bases) {
        if let Some(solver) = slot {
            *solver.base_mut() = base;
        }
    }
}

impl SessionObjects<'_> {
    /// `m_apUserPointer[0..2]` of every rail vertex use's intersection point: the
    /// surfaces of the faces of its geometry's solver's edge uses 0 and 1, as
    /// `SmFilletGeom::MakeRailEdgeuseOn` sets them.
    fn tsect_surfaces(&self) -> Result<BTreeMap<VertexuseId, [Option<Surface>; 2]>> {
        let mut out = BTreeMap::new();
        for (vu, _) in self.pseudo.fillet_vertexuses() {
            let VertexuseOwner::Edgeuse(Some(eu)) = self.pseudo.brep().vertexuse(vu)?.owner
            else {
                continue;
            };
            let Some(edge) = self.pseudo.brep().edgeuse(eu)?.edge else {
                continue;
            };
            let is_rail = self
                .pseudo
                .fillet_edge(edge)
                .is_ok_and(|e| e.fillet_edge_type() == FilletEdgeType::Rail);
            if !is_rail {
                continue;
            }
            let Some(geom) = self
                .pseudo
                .fillet_edgeuse(eu)
                .ok()
                .and_then(crate::blend::fillet_geom::p1::FilletEdgeuse::fillet_geom)
            else {
                continue;
            };
            let Some(solver) = self.pseudo.geom(geom)?.solver else {
                continue;
            };
            let Some(Some(solver)) = self.solvers.get(solver.0 as usize) else {
                continue;
            };
            let mut surfaces = [None, None];
            for (index, slot) in surfaces.iter_mut().enumerate() {
                if let Some(base_eu) = solver.base().edgeuse(index)? {
                    let face = crate::topology::edgeuse::p1::face(&*self.target, base_eu)?;
                    if let Some(face) = face {
                        // (the face's own surface object: its owner is kept)
                        *slot = self
                            .target
                            .face(face)?
                            .surface
                            .as_ref()
                            .map(Surface::clone_keeping_owner);
                    }
                }
            }
            out.insert(vu, surfaces);
        }
        Ok(out)
    }

    /// Moves the executive's members, the solvers' bases and the corner table into
    /// corner form, `id` taken out of the table. Nothing moves when a conversion fails.
    fn take_corner_state(
        &mut self,
        executive: &FilletExecutive<FilletObject>,
        id: FilletCornerId,
        detach_solvers: bool,
    ) -> Result<Taken> {
        // conversions first: an unconvertible member leaves the session untouched
        let mut edgeuse_solver_map = BTreeMap::new();
        for (&eu, &solver) in &executive.edgeuse_solver_map {
            let (Some(eu), FilletObject::Solver(solver)) = (target_edgeuse(eu), solver) else {
                return Err(invalid("edge-use solver map entry is not a target edge use"));
            };
            edgeuse_solver_map.insert(eu, solver);
        }
        let mut topo_edges = Vec::with_capacity(executive.topo_edges.len());
        for entry in &executive.topo_edges {
            topo_edges.push(match entry {
                None => None,
                Some(object) => Some(
                    pseudo_edge(*object)
                        .ok_or_else(|| invalid("topology edge is not a pseudo-B-rep edge"))?,
                ),
            });
        }
        let mut vertex_corner_map = BTreeMap::new();
        for (vertex, entry) in &executive.vertex_corner_map {
            if let (Some(v), VertexEntry::Corner(FilletObject::Corner(c))) =
                (target_vertex(*vertex), entry)
            {
                vertex_corner_map.insert(v, *c);
            }
        }
        let fillet_corners = executive
            .corners
            .iter()
            .map(|c| match c {
                FilletObject::Corner(c) => Ok(*c),
                _ => Err(invalid("executive corner is not a session corner")),
            })
            .collect::<Result<Vec<_>>>()?;
        let fillet_solvers = executive
            .solvers
            .iter()
            .map(|s| match s {
                FilletObject::Solver(s) => Ok(*s),
                _ => Err(invalid("executive solver is not a session solver")),
            })
            .collect::<Result<Vec<_>>>()?;
        let index = id.0 as usize;
        if index >= self.corners.len() || u32::try_from(self.corners.len()).is_err() {
            return Err(invalid("corner is not live"));
        }

        // then the moves, which cannot fail
        let mut corners: BTreeMap<FilletCornerId, FilletCorner> = BTreeMap::new();
        let mut running = None;
        for (i, corner) in std::mem::take(&mut self.corners).into_iter().enumerate() {
            let cid = FilletCornerId(u32::try_from(i).unwrap_or(u32::MAX));
            if cid == id {
                running = Some(corner);
            } else {
                corners.insert(cid, corner);
            }
        }
        let solvers = if detach_solvers {
            detach_bases(&mut self.solvers)
        } else {
            Vec::new()
        };
        let data = FilletExecutiveData {
            edgeuse_solver_map,
            extended_surfaces_map: std::mem::take(&mut self.extended_surfaces),
            do_global_merge: executive.options.global_merge,
            topo_edges,
            error: ErrorInfo::default(),
            corners,
            fillet_corners,
            vertex_corner_map,
            fillet_solvers,
            ..FilletExecutiveData::default()
        };
        Ok(Taken {
            // (index checked above, so the running corner was found)
            corner: running.ok_or_else(|| invalid("corner is not live"))?,
            data,
            solvers,
        })
    }

    /// Moves everything [`SessionObjects::take_corner_state`] took back, with the
    /// corner code's changes. Infallible in its moves; the members the session keeps
    /// in its own form are converted back.
    fn restore_corner_state(
        &mut self,
        executive: &mut FilletExecutive<FilletObject>,
        id: FilletCornerId,
        taken: Taken,
    ) -> Result<()> {
        let Taken {
            corner,
            mut data,
            solvers,
        } = taken;
        attach_bases(&mut self.solvers, solvers);
        data.corners.insert(id, corner);
        let expected = data.corners.len();
        self.corners = std::mem::take(&mut data.corners).into_values().collect();
        self.extended_surfaces = std::mem::take(&mut data.extended_surfaces_map);
        executive.options.global_merge = data.do_global_merge;
        executive.topo_edges = data
            .topo_edges
            .iter()
            .map(|e| e.map(|e| FilletObject::Topology(FilletBody::Pseudo, EntityRef::Edge(e))))
            .collect();
        if data.error != ErrorInfo::default() {
            executive.error = error_info(&data.error);
        }
        // the corner table is indexed by id: new or removed corners must keep it dense
        let dense = self.corners.len() == expected
            && data
                .fillet_corners
                .iter()
                .all(|c| (c.0 as usize) < self.corners.len());
        executive.corners = data
            .fillet_corners
            .iter()
            .map(|c| FilletObject::Corner(*c))
            .collect();
        executive.solvers = data
            .fillet_solvers
            .iter()
            .map(|s| FilletObject::Solver(*s))
            .collect();
        for (v, c) in &data.vertex_corner_map {
            executive.vertex_corner_map.insert(
                Some(FilletObject::Topology(FilletBody::Target, EntityRef::Vertex(*v))),
                VertexEntry::Corner(FilletObject::Corner(*c)),
            );
        }
        if !data.deleted_geoms.is_empty() {
            return Err(invalid(
                "a corner call deleted fillet geometries the session does not track",
            ));
        }
        if dense {
            Ok(())
        } else {
            Err(invalid("a corner call left the corner table sparse"))
        }
    }

    /// Runs `run` on corner `id`, taken out of the corner table, with the executive data
    /// the corner lookups read ([`FilletExecutiveData`]: the other corners by id, the
    /// extended-surface map, the maps, the global-merge switch) and the session's other
    /// objects, the solvers untouched; then moves everything back, also when `run` fails.
    ///
    /// # Errors
    ///
    /// As [`SessionObjects::with_corner_env`].
    pub(super) fn with_corner_data<T>(
        &mut self,
        executive: &mut FilletExecutive<FilletObject>,
        id: FilletCornerId,
        run: impl FnOnce(&mut Self, &mut FilletCorner, &mut FilletExecutiveData) -> Result<T>,
    ) -> Result<T> {
        let mut taken = self.take_corner_state(executive, id, false)?;
        let result = run(self, &mut taken.corner, &mut taken.data);
        let restored = self.restore_corner_state(executive, id, taken);
        let value = result?;
        restored?;
        Ok(value)
    }

    /// Runs `run` on corner `id` with a [`CornerEnv`] over this session's objects
    /// (no fillet B-rep: `m_pFilletBrep` is `NULL` before `CreateFilletBrep`), then
    /// moves every member back, also when `run` fails.
    ///
    /// # Errors
    ///
    /// A session member the corner form cannot hold (see the module docs); those of
    /// `run`; a corner call that changed what the session cannot take back.
    pub(super) fn with_corner_env<T>(
        &mut self,
        k: &mut Kernel,
        executive: &mut FilletExecutive<FilletObject>,
        id: FilletCornerId,
        run: impl FnOnce(&mut FilletCorner, &mut CornerEnv<'_>) -> Result<T>,
    ) -> Result<T> {
        let tsect_surfaces = self.tsect_surfaces()?;
        let mut taken = self.take_corner_state(executive, id, true)?;
        let result = {
            let Taken {
                corner,
                data,
                solvers,
            } = &mut taken;
            let mut env = CornerEnv {
                k,
                target: &mut *self.target,
                pseudo: &mut self.pseudo,
                exec: data,
                corner: Some(id),
                tsect_surfaces: &tsect_surfaces,
                solvers: solvers.as_mut_slice(),
                fillet: None,
            };
            run(corner, &mut env)
        };
        let restored = self.restore_corner_state(executive, id, taken);
        let value = result?;
        restored?;
        Ok(value)
    }
}
