//! Topology at the end of one filleted edge.
//!
//! Upstream: `source/SMLib/src/SmFilletCorner.cpp` (44–119, 443–506, 5990–6443).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`FilletCorner::make_nx1_corner_topology`], [`FilletCorner::make_nx1_corner_topology_with`] | `SmFilletNx1Corner::MakeCornerTopology` |
//! | `find_side_edgeuse`, `find_closed_side_edgeuse` | `sm_FindSideEdgeuse`, `sm_FindClosedSideEdgeuse` |
//! | `same_convexity_edges` | `same_convexity_edges` |
//! | `find_surface_extension_dist_at_corner`, [`CornerExtensionRadius`] | `sm_FindSurfaceExtensionDistAtCorner` |
//! | [`Nx1MixedEnv`] | the executive's `m_vExtendedSurfacesMap`, `AddTopoEdges` |
//!
//! Deviations: the caller lends the actual target, pseudo-B-rep and solver, and for
//! a mixed-convexity corner the solver's radius data and the executive's
//! extended-surface map and topological-edge list ([`Nx1MixedEnv`]); without them
//! that corner stops (`Unported`) after its native vertex/rail prefix, as before.
//! Missing topology that native dereferences is an input error. A variable radius
//! law that cannot be evaluated fails the corner (upstream reads its status as the
//! radius). The multi-neighbour end face (`bMultiNeighborEndFace`) stops `Unported`
//! after its topology, at the corner geometry and coincidence check. There is no
//! recursion.

use std::collections::BTreeMap;

use crate::blend::fillet_corner::p1::face_by_two_edgeuses;
use crate::blend::fillet_corner::p3::{FilletCorner, FilletCornerKind};
use crate::blend::fillet_geom::p1::{
    FilletBrep, FilletCornerId, FilletEdgeType, FilletVertex, FilletVertexType,
};
use crate::blend::fillet_solver::p1::FilletSolver;
use crate::error::{EntityRef, Result, SolidError};
use crate::geom::geom_utility::p2::ContinuityType;
use crate::intersect::point_class::PointClassType;
use crate::kernel::Kernel;
use crate::math::{Point3, EFF_ZERO};
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::surface::surface::Surface;
use crate::topology::edge;
use crate::topology::edgeuse::{p1 as eu1, p2 as eu2};
use crate::topology::{Brep, EdgeId, EdgeuseId, FaceId, Vertex, VertexId};

const OP: &str = "blend::fillet_corner::topology::make_nx1_corner_topology";

/// Actual objects used by `SmFilletNx1Corner::MakeCornerTopology`.
#[derive(Debug)]
pub struct Nx1TopologyEnv<'a> {
    /// The executive's kernel.
    pub k: &'a Kernel,
    /// The target whose edge uses the corner stores.
    pub target: &'a mut Brep,
    /// The executive's actual pseudo-B-rep.
    pub pseudo: &'a mut FilletBrep,
    /// This corner's identity in the executive's corner store.
    pub corner: FilletCornerId,
    /// The loaded solver of the corner's filleted edge use.
    pub solver: &'a FilletSolver,
    /// The executive's live global-merge switch.
    pub do_global_merge: &'a mut bool,
}

/// The solver data `sm_FindSurfaceExtensionDistAtCorner` reads (by the solver's
/// type): the derived solver's radii, distance or radius law.
pub enum CornerExtensionRadius<'a> {
    /// `SM_FS_CONST_RADIUS` / `SM_FS_CONST_RADIUS_ASSISTED`: `GetFilletRadius(0)`,
    /// `GetFilletRadius(1)`.
    ConstRadius([f64; 2]),
    /// `SM_FS_CONST_DIST`: `GetDistance()`.
    ConstDist(f64),
    /// `SM_FS_VARIABLE_RADIUS`: `GetFilletRadius(dParameter)`, the law's radius at an
    /// edge parameter.
    VariableRadius(&'a mut dyn FnMut(f64) -> Result<f64>),
    /// Any other solver type (`SE(SM_ERR)`: logged, the distance stays 0).
    Other,
}

impl core::fmt::Debug for CornerExtensionRadius<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::ConstRadius(r) => f.debug_tuple("ConstRadius").field(r).finish(),
            Self::ConstDist(d) => f.debug_tuple("ConstDist").field(d).finish(),
            Self::VariableRadius(_) => f.write_str("VariableRadius(..)"),
            Self::Other => f.write_str("Other"),
        }
    }
}

/// What the mixed-convexity branch of `SmFilletNx1Corner::MakeCornerTopology` reads
/// and writes beyond [`Nx1TopologyEnv`]: the solver's radius data, the executive's
/// `m_vExtendedSurfacesMap` (`AddExtendedOriginalSurface`) and its `m_vTopoEdges`
/// (`AddTopoEdges`, the pseudo-B-rep fillet edges appended in order).
#[derive(Debug)]
pub struct Nx1MixedEnv<'a> {
    /// The solver's radius data for the extension distance.
    pub radius: CornerExtensionRadius<'a>,
    /// `m_vExtendedSurfacesMap`.
    pub extended_surfaces: &'a mut BTreeMap<FaceId, Surface>,
    /// The fillet edges `AddTopoEdges` registers, in order.
    pub topo_edges: &'a mut Vec<EdgeId>,
}

fn invalid(reason: &'static str) -> SolidError {
    SolidError::invalid_input(OP, reason)
}

/// `sm_FindSurfaceExtensionDistAtCorner(pFilEdgeuse, pFilSolver, cpNx1CornerVert)`: an
/// extension distance for the surface at an Nx1 corner, 7.5 times the largest constant
/// radius, the constant distance or the variable radius at the corner's end of the
/// filleted edge (0 for the other solver types).
///
/// # Errors
/// The variable radius law's evaluation (upstream returns its status as the radius);
/// the topology queries.
fn find_surface_extension_dist_at_corner(
    target: &Brep,
    filleted: EdgeuseId,
    radius: &mut CornerExtensionRadius<'_>,
    corner_vertex: VertexId,
) -> Result<f64> {
    let ext_dist = match radius {
        // return largest fillet radius
        CornerExtensionRadius::ConstRadius(radii) => {
            let mut ext_dist = 0.0;
            for &rad in radii.iter() {
                if rad > ext_dist {
                    ext_dist = rad;
                }
            }
            ext_dist
        }
        // return constant distance
        CornerExtensionRadius::ConstDist(distance) => *distance,
        // return radius at vertex
        CornerExtensionRadius::VariableRadius(radius_at) => {
            let ivl = target.edge(edge_of(target, filleted)?)?.interval;
            let vert = vertex_of(target, filleted)?;
            let orientation = eu1::orientation(target, filleted)?;
            let at_max = (Some(corner_vertex) == vert && orientation == OrientType::Opposite)
                || (Some(corner_vertex) != vert && orientation == OrientType::Same);
            radius_at(if at_max { ivl.max() } else { ivl.min() })?
        }
        CornerExtensionRadius::Other => 0.0,
    };

    // In prog_test, Fillet test my_test_suite_0 iter: 2 will create a surface with a
    // slight crease for extension values greater than about 7.5.
    let rad_factor = 7.5;
    Ok(rad_factor * ext_dist)
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

fn find_closed_side_edgeuse(
    target: &Brep,
    solver: &FilletSolver,
    side: EdgeuseId,
    rail: usize,
) -> Result<EdgeuseId> {
    let filleted = need(solver.edgeuse(0)?, "solver without first edge use")?;
    let same_vertex = vertex_of(target, side)? == vertex_of(target, filleted)?;
    let neighbour = if same_vertex {
        need(
            eu1::cw_edgeuse(target, filleted)?,
            "filleted edge use without CW neighbour",
        )?
    } else {
        need(
            eu1::ccw_edgeuse(target, filleted)?,
            "filleted edge use without CCW neighbour",
        )?
    };
    if (rail == 0) == same_vertex {
        need(eu1::mate(target, neighbour)?, "side edge use without mate")
    } else {
        Ok(neighbour)
    }
}

fn same_convexity_edges(
    k: &Kernel,
    target: &mut Brep,
    first: EdgeuseId,
    second: EdgeuseId,
    vertex: VertexId,
) -> Result<bool> {
    let concave = if eu2::is_convex_radial_sector(k, target, first, 5, false, 0.0)? {
        if edge::is_manifold(target, edge_of(target, first)?)? {
            let mate = need(eu1::mate(target, first)?, "filleted edge use without mate")?;
            Some(radial(target, mate)?)
        } else {
            None
        }
    } else {
        Some(first)
    };
    // Native SER(SM_ERR) returns a nonzero status through this SmBoolean.
    let Some(concave) = concave else {
        return Ok(true);
    };
    let second_edge = eu1::edge(target, second)?;
    let mut side = need(
        eu1::corner_mate_edgeuse(target, concave, Some(vertex))?,
        "concave edge use without corner mate",
    )?;
    if eu1::edge(target, side)? != second_edge {
        side = need(
            eu1::corner_mate_edgeuse(target, radial(target, concave)?, Some(vertex))?,
            "concave radial without corner mate",
        )?;
    }
    Ok(!eu2::is_convex_radial_sector(
        k, target, side, 5, false, 0.0,
    )?)
}

fn is_fillet_side_face(target: &Brep, filleted: EdgeuseId, face: Option<FaceId>) -> Result<bool> {
    if face == eu1::face(target, filleted)? {
        return Ok(true);
    }
    Ok(face == eu1::face(target, radial(target, filleted)?)?)
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
    /// Construct the rails and end rib of an actual `NByOne` corner.
    ///
    /// Allocations and corner-list insertions occur in native order and remain
    /// visible if a later operation fails. Vertex geometry is computed later.
    ///
    /// # Errors
    /// Invalid topology or a different corner kind; those of the actual solver,
    /// convexity and pseudo-topology operations. A mixed-convexity corner is
    /// `Unported` after creating its native rail/vertex prefix: it needs the
    /// executive's data, [`FilletCorner::make_nx1_corner_topology_with`].
    pub fn make_nx1_corner_topology(&mut self, env: Nx1TopologyEnv<'_>) -> Result<()> {
        self.make_nx1_corner_topology_with(env, None)
    }

    /// [`FilletCorner::make_nx1_corner_topology`] with the solver's radius data and the
    /// executive's extended-surface map and topological-edge list a mixed-convexity
    /// corner needs (`mixed`): the end rib, the end face's extension and, at the
    /// corner's own vertex, a fillet vertex and the two fillet edges to the rails' ends
    /// (`SmFilletNx1Corner::MakeCornerTopology`, SmFilletCorner.cpp:6150-6442).
    ///
    /// # Errors
    /// Those of [`FilletCorner::make_nx1_corner_topology`]; at a mixed-convexity corner
    /// without `mixed`, `Unported`; those of the extension distance, the surface
    /// extension and copy; an end face between more than one neighbour face
    /// (`bMultiNeighborEndFace`) is `Unported` after its topology, at the corner
    /// geometry and coincidence check upstream runs there.
    #[allow(clippy::too_many_lines)] // one upstream routine
    pub fn make_nx1_corner_topology_with(
        &mut self,
        env: Nx1TopologyEnv<'_>,
        mut mixed_env: Option<Nx1MixedEnv<'_>>,
    ) -> Result<()> {
        let Nx1TopologyEnv {
            k,
            target,
            pseudo,
            corner,
            solver,
            do_global_merge,
        } = env;
        if self.kind != FilletCornerKind::NByOne {
            return Err(invalid("requires an NByOne corner"));
        }
        let filleted = need(
            self.solver_eus
                .as_deref()
                .and_then(|eus| eus.first())
                .copied(),
            "corner without filleted edge use",
        )?;
        let geometry = solver.first_fillet_geom()?;
        let all = need(self.all_eus.as_deref(), "corner without outward edge uses")?;
        let total = all.len();
        let index = need(
            all.iter().position(|&eu| eu == filleted),
            "filleted edge use not in corner",
        )?;
        let mut sides = [
            *need(
                all.get((index + total - 1) % total),
                "preceding side edge use",
            )?,
            *need(all.get((index + 1) % total), "following side edge use")?,
        ];
        let [side0, side1] = sides;
        let closed_side = side0 == side1;
        let vertex = need(self.vertex, "corner without vertex")?;
        let same = [
            same_convexity_edges(k, target, filleted, side0, vertex)?,
            same_convexity_edges(k, target, filleted, side1, vertex)?,
        ];
        let [same0, same1] = same;
        let mixed = !same0 || !same1;
        // when working with mixed convexity: the side edge to be extended, and the
        // face to extend when the end face's surface can be a surface extension
        let mut extend_edge: Option<EdgeId> = None;
        let mut end_face = None;

        for rail in 0..2 {
            let side = *need(sides.get(rail), "side index")?;
            let fv = new_corner_vertex(k, pseudo, corner)?;
            self.vertices.push(fv);
            let (side_use, solver_rail) = if closed_side {
                (find_closed_side_edgeuse(target, solver, side, rail)?, rail)
            } else {
                let found = need(
                    find_side_edgeuse(target, filleted, side)?,
                    "no side edge use intersects rail",
                )?;
                (
                    found,
                    solver.find_index_of_rail_x_side_edgeuse(target, found)?,
                )
            };
            pseudo
                .fillet_vertex_mut(fv)?
                .set_point_class(PointClassType::Edgeuse, Some(EntityRef::Edgeuse(side_use)));
            pseudo.make_rail_edgeuse_on(
                geometry,
                solver_rail,
                Some(fv),
                OrientType::Unknown,
                None,
                crate::blend::fillet_geom::p1::CornerLookup {
                    this: Some((Some(corner), self)),
                    ..Default::default()
                },
                solver,
                target,
            )?;

            let face = eu1::face(target, side)?;
            if mixed && *need(same.get(rail), "convexity index")? {
                pseudo
                    .fillet_vertex_mut(fv)?
                    .set_fillet_vertex_type(FilletVertexType::RailXExtendedEdgeuse);
                // extend the edge with the same convexity: it shares a face with the
                // rail-to-be it will intersect; the end face is the one between the
                // extended edge and the other side edge
                let other = *need(sides.get(1 - rail), "other side index")?;
                end_face = face_by_two_edgeuses(target, Some(side), Some(other), None)?;
                extend_edge = Some(edge_of(target, side)?);
                if !is_fillet_side_face(target, filleted, face)? {
                    *need(sides.get_mut(rail), "side index")? = radial(target, side)?;
                }
            } else {
                pseudo
                    .fillet_vertex_mut(fv)?
                    .set_fillet_vertex_type(FilletVertexType::RailXEdgeuse);
                if is_fillet_side_face(target, filleted, face)? {
                    *need(sides.get_mut(rail), "side index")? = radial(target, side)?;
                }
            }
        }
        if !mixed && total > 3 {
            *do_global_merge = true;
            return Ok(());
        }

        // compute how much to extend the end surface to make sure it covers the 'hole'
        let mut ext_dist = if mixed {
            let Some(mixed_env) = mixed_env.as_mut() else {
                // UNPORTED(blend.fillet_corner_p1): sm_FindSurfaceExtensionDistAtCorner
                // without the executive's data (make_nx1_corner_topology_with)
                return Err(SolidError::unported(
                    "sm_FindSurfaceExtensionDistAtCorner",
                    "blend.fillet_corner_p1",
                ));
            };
            find_surface_extension_dist_at_corner(target, filleted, &mut mixed_env.radius, vertex)?
        } else {
            0.0
        };

        let start = self.vertices.first().copied();
        let end = self.vertices.get(1).copied();
        let rib = pseudo.make_fillet_edge(start, end, Some(corner))?;
        pseudo
            .fillet_edge_mut(rib)?
            .set_fillet_edge_type(FilletEdgeType::FilletXSideFace);
        let primary = need(
            eu1::primary_edgeuse(pseudo.brep(), rib)?,
            "end rib without primary edge use",
        )?;
        pseudo.set_edgeuse_fillet_geom(primary, Some(geometry))?;
        self.edges.push(rib);

        // In the mixed-convexity case the end face is still unknown when both side
        // edges have opposite convexity to the filleted edge: the face connecting the
        // two side edges is used (or extended).
        let mut multi_neighbor_end_face = false;
        if end_face.is_none() {
            let [first, second] = sides;
            end_face = face_by_two_edgeuses(target, Some(first), Some(second), None)?;
            if end_face.is_none() {
                // more than 3 edges at the corner, the neighbour edges connected by a
                // sequence of faces (only at a mixed-convexity corner here)
                if !mixed || total <= 3 {
                    return Err(SolidError::sm_err(OP));
                }
                multi_neighbor_end_face = true;
                // the cliff edge: the side edge with the filleted edge's convexity,
                // else (both cliff edges) the first
                let pick = usize::from(!same0 && same1);
                let cliff = *need(sides.get(pick), "cliff side index")?;
                end_face = eu1::face(target, cliff)?;
                extend_edge = Some(edge_of(target, cliff)?);
            }
        }
        pseudo.fillet_edge_mut(rib)?.set_original_face(end_face);

        // compute how much to extend the side surface to make sure it covers the 'hole'
        let mut do_extend = extend_edge.is_some();
        if ext_dist > EFF_ZERO
            && !do_extend
            && mixed_env
                .as_ref()
                .is_some_and(|m| matches!(m.radius, CornerExtensionRadius::ConstRadius(_)))
        {
            ext_dist *= 0.1; // use only a small extension
            do_extend = true; // for constant radius only
        }

        // When a surface extension is needed to define the end face's shape, build and
        // save it.
        let mut new_extend_surface_kept = false;
        if do_extend {
            let mixed_env = mixed_env
                .as_mut()
                .ok_or_else(|| invalid("surface extension without the executive's map"))?;
            let face = need(end_face, "extension without end face")?;
            let surface = target
                .face(face)?
                .surface
                .as_ref()
                .ok_or_else(|| SolidError::null_pointer(OP))?;
            // create an extended surface large enough to fill in the end cap
            let new_surface = surface.create_extended_surface(k, ext_dist, ContinuityType::G1)?;
            // Copy it, when possible as an analytic surface
            let new_extend_surface = new_surface.copy_and_add_analytics(k)?;
            // register the face / extended surface pair; an existing larger one stays
            new_extend_surface_kept = match mixed_env.extended_surfaces.get(&face) {
                Some(old) => !new_extend_surface
                    .natural_uv_domain()
                    .is_contained_by(&old.natural_uv_domain(), EFF_ZERO),
                None => true,
            };
            FilletCorner::add_extended_original_surface_with_map(
                &mut *mixed_env.extended_surfaces,
                face,
                new_extend_surface,
            );
        }

        // When adding an end surface to fill the 'hole', set up the end cap's topology.
        if mixed {
            let mixed_env = need(mixed_env.as_mut(), "mixed corner without the executive")?;
            // add a fillet vertex at the original vertex's location
            let fv = new_corner_vertex(k, pseudo, corner)?;
            pseudo
                .fillet_vertex_mut(fv)?
                .set_point_class(PointClassType::Vertex, Some(EntityRef::Vertex(vertex)));
            pseudo
                .fillet_vertex_mut(fv)?
                .set_fillet_vertex_type(FilletVertexType::OnVertex);
            self.vertices.push(fv);

            // the fillet edges from the expected rail 2 / side edge 2 intersection to
            // the original vertex, then from it to the rail 1 / side edge intersection
            for (from, to, side_index) in [(1, 2, 1), (2, 0, 0)] {
                let start = self.vertices.get(from).copied();
                let end = self.vertices.get(to).copied();
                let new_edge = pseudo.make_fillet_edge(start, end, Some(corner))?;
                let side = *need(sides.get(side_index), "side index")?;
                // the edge type follows the analysis saved in the rail's fillet vertex
                let rail_vertex = need(self.vertices.get(side_index).copied(), "rail vertex")?;
                let extended = pseudo.fillet_vertex(rail_vertex)?.fillet_vertex_type()
                    == FilletVertexType::RailXExtendedEdgeuse;
                if extended {
                    let fe = pseudo.fillet_edge_mut(new_edge)?;
                    fe.set_fillet_edge_type(FilletEdgeType::OnExtendEdge);
                    fe.set_original_edge(extend_edge);
                } else {
                    pseudo
                        .fillet_edge_mut(new_edge)?
                        .set_fillet_edge_type(FilletEdgeType::OnEdge);
                    if total <= 3 {
                        // registered, processed when removing topological vertices
                        mixed_env.topo_edges.push(new_edge);
                    }
                    let side_edge = edge_of(target, side)?;
                    pseudo
                        .fillet_edge_mut(new_edge)?
                        .set_original_edge(Some(side_edge));
                }
                let side_face = eu1::face(target, side)?;
                pseudo
                    .fillet_edge_mut(new_edge)?
                    .set_original_face(side_face);
                self.edges.push(new_edge);
            }

            // all the corner topology is in place; when the end face connects to more
            // than the fillet surface and one other face, upstream sets the corner's
            // vertex positions and side edge curves here and checks the extended
            // surface against the other side edge
            if multi_neighbor_end_face {
                let _ = new_extend_surface_kept;
                // UNPORTED(blend.fillet_corner_p1): CalcCornerVertGeom, CalcCornerEdgeGeom
                // and SimpleCoincidenceChecker of the multi-neighbour end face
                return Err(SolidError::unported(
                    "SmFilletNx1Corner::MakeCornerTopology (multi-neighbour end face)",
                    "blend.fillet_corner_p1",
                ));
            }
        }
        Ok(())
    }
}
