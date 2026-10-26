//! Topology of the 3x2 mixed-convexity fillet corner.
//!
//! Upstream: `source/SMLib/src/SmFilletCorner.cpp` (443–506, 4755–4897, 5990 ff.).
//!
//! | Rust | Upstream |
//! |------|----------|
//! | [`FilletCorner::make_3x2_mixed_corner_topology`] | `SmFillet3x2MixedConvexityCorner::MakeCornerTopology` |
//! | [`Mixed3x2TopologyEnv`] | the executive's pseudo-B-rep, `GetFilletSolverOfEdgeuse` and `m_vExtendedSurfacesMap` |
//! | `find_side_edgeuse` (private) | `sm_FindSideEdgeuse` (as `topology.rs`'s) |
//! | `find_surface_extension_dist_at_corner` (private) | `sm_FindSurfaceExtensionDistAtCorner` (as `topology.rs`'s) |
//!
//! Deviations: the caller lends the target, the pseudo-B-rep, the solvers of the two
//! filleted edge uses, the convex solver's radius data and the executive's
//! extended-surface map. Missing topology that upstream dereferences (no side edge use,
//! no common face) is an input error. A variable radius law that cannot be evaluated
//! fails the corner (upstream reads its status as the radius). There is no recursion.

use std::collections::BTreeMap;

use crate::blend::fillet_corner::p1::face_by_two_edgeuses;
use crate::blend::fillet_corner::p3::{FilletCorner, FilletCornerKind};
use crate::blend::fillet_corner::topology::CornerExtensionRadius;
use crate::blend::fillet_geom::p1::{
    CornerLookup, FilletBrep, FilletCornerId, FilletEdgeType, FilletVertex, FilletVertexType,
};
use crate::blend::fillet_solver::p1::FilletSolver;
use crate::error::{EntityRef, Result, SolidError};
use crate::geom::geom_utility::p2::ContinuityType;
use crate::intersect::point_class::PointClassType;
use crate::kernel::Kernel;
use crate::math::Point3;
use crate::nurbs::b_spline_surface_n_lib::p2::OrientType;
use crate::surface::surface::Surface;
use crate::topology::edgeuse::{p1 as eu1, p2 as eu2};
use crate::topology::{Brep, EdgeId, EdgeuseId, FaceId, Vertex, VertexId};

const OP: &str = "blend::fillet_corner::mixed_3x2::make_3x2_mixed_corner_topology";

/// What `SmFillet3x2MixedConvexityCorner::MakeCornerTopology` works on.
#[derive(Debug)]
pub struct Mixed3x2TopologyEnv<'a> {
    /// The executive's kernel.
    pub k: &'a Kernel,
    /// The B-rep being filleted.
    pub target: &'a mut Brep,
    /// The executive's pseudo-B-rep.
    pub pseudo: &'a mut FilletBrep,
    /// This corner's identity in the executive's corner store.
    pub corner: FilletCornerId,
    /// The loaded solvers of the two filleted edge uses, `(*m_pSolverEUs)[0]` and `[1]`
    /// (`GetFilletSolverOfEdgeuse`).
    pub solvers: [&'a FilletSolver; 2],
    /// The radius data of each of those solvers, for the extension distance (the convex
    /// one's is read).
    pub radii: [CornerExtensionRadius<'a>; 2],
    /// The executive's `m_vExtendedSurfacesMap` (`AddExtendedOriginalSurface`).
    pub extended_surfaces: &'a mut BTreeMap<FaceId, Surface>,
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

fn radial(target: &Brep, eu: EdgeuseId) -> Result<EdgeuseId> {
    need(eu1::radial(target, eu)?, "edge use without radial")
}

fn vertex_of(target: &Brep, eu: EdgeuseId) -> Result<Option<VertexId>> {
    let vu = need(eu1::vertexuse(target, eu)?, "edge use without vertex use")?;
    Ok(target.vertexuse(vu)?.vertex)
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
    let side = radial(target, side)?;
    let ccw = need(
        eu1::ccw_edgeuse(target, side)?,
        "side radial without CCW neighbour",
    )?;
    Ok((filleted_edge == eu1::edge(target, ccw)?).then_some(side))
}

/// `sm_FindSurfaceExtensionDistAtCorner(pFilEdgeuse, pFilSolver, cpCornerVert)`: 7.5
/// times the largest constant radius, the constant distance or the variable radius at
/// the corner's end of the filleted edge (0 for the other solver types).
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
    /// `SmFillet3x2MixedConvexityCorner::MakeCornerTopology`: three edges at the
    /// corner, a convex and a concave one filleted. The concave fillet ends in a cross
    /// section from its rail on the side edge to where both fillets' rails meet in their
    /// common face; the convex fillet's cross section runs from there to its rail on the
    /// (extended) side edge; two edges along the side edge, through a vertex on the
    /// corner's own vertex, close the end. The surface of the face between the side edge
    /// and the concave edge is extended (by the convex solver's extension distance) and
    /// registered as that face's extended original surface.
    ///
    /// # Errors
    ///
    /// Invalid topology or a different corner kind; `SM_ERR` unless the corner has three
    /// edge uses; those of the convexity test, the extension distance, the surface
    /// extension and copy, and the pseudo-topology operations.
    #[allow(clippy::too_many_lines)] // one upstream routine
    pub fn make_3x2_mixed_corner_topology(&mut self, env: Mixed3x2TopologyEnv<'_>) -> Result<()> {
        let Mixed3x2TopologyEnv {
            k,
            target,
            pseudo,
            corner,
            solvers,
            radii,
            extended_surfaces,
        } = env;
        if self.kind != FilletCornerKind::ThreeByTwoMixed {
            return Err(invalid("requires a 3x2 mixed-convexity corner"));
        }
        let all = need(self.all_eus.clone(), "corner without outward edge uses")?;
        if all.len() != 3 {
            return Err(SolidError::sm_err(OP));
        }
        let solver_eu = |i: usize| {
            need(
                self.solver_eus
                    .as_deref()
                    .and_then(|eus| eus.get(i))
                    .copied(),
                "corner without its filleted edge uses",
            )
        };
        let (mut convex, mut concave) = (solver_eu(0)?, solver_eu(1)?);
        let [mut convex_solver, mut concave_solver] = solvers;
        let [mut convex_radius, mut concave_radius] = radii;
        if !eu2::is_convex_radial_sector(k, target, concave, 5, false, 0.0)? {
            core::mem::swap(&mut concave, &mut convex);
            core::mem::swap(&mut concave_solver, &mut convex_solver);
            core::mem::swap(&mut concave_radius, &mut convex_radius);
        }
        let convex_fg = convex_solver.first_fillet_geom()?;
        let concave_fg = concave_solver.first_fillet_geom()?;
        let vertex = need(self.vertex, "corner without vertex")?;

        let curr = need(
            all.iter().copied().find(|&eu| eu != convex && eu != concave),
            "no side edge use",
        )?;

        // This edge use is 'in-between' two filleted edge uses
        let fv0 = new_corner_vertex(k, pseudo, corner)?;
        let side = need(
            find_side_edgeuse(target, concave, curr)?,
            "no side edge use intersects the concave rail",
        )?;
        let rail = concave_solver.find_index_of_rail_x_side_edgeuse(target, side)?;
        pseudo.make_rail_edgeuse_on(
            concave_fg,
            rail,
            Some(fv0),
            OrientType::Unknown,
            None,
            lookup(corner, self),
            concave_solver,
            target,
        )?;
        let record = pseudo.fillet_vertex_mut(fv0)?;
        record.set_point_class(PointClassType::Edgeuse, Some(EntityRef::Edgeuse(side)));
        record.set_fillet_vertex_type(FilletVertexType::RailXEdgeuse);
        self.vertices.push(fv0);

        let common_face = face_by_two_edgeuses(target, Some(convex), Some(concave), Some(curr))?;
        let fv1 = new_corner_vertex(k, pseudo, corner)?;
        let record = pseudo.fillet_vertex_mut(fv1)?;
        record.set_point_class(PointClassType::Face, common_face.map(EntityRef::Face));
        record.set_fillet_vertex_type(FilletVertexType::RailXRail);
        self.vertices.push(fv1);

        // Setup rails of fillet
        for (fg, rail_solver) in [(concave_fg, concave_solver), (convex_fg, convex_solver)] {
            let rail = rail_solver.find_index_of_rail_on_face(target, common_face)?;
            pseudo.make_rail_edgeuse_on(
                fg,
                rail,
                Some(fv1),
                OrientType::Unknown,
                None,
                lookup(corner, self),
                rail_solver,
                target,
            )?;
        }

        // Make a cross section on the concave fillet
        let fe0 = pseudo.make_fillet_edge(Some(fv0), Some(fv1), Some(corner))?;
        pseudo
            .fillet_edge_mut(fe0)?
            .set_fillet_edge_type(FilletEdgeType::CrossSection);
        let primary = need(
            eu1::primary_edgeuse(pseudo.brep(), fe0)?,
            "fillet edge without primary edge use",
        )?;
        pseudo.set_edgeuse_fillet_geom(primary, Some(concave_fg))?;
        self.edges.push(fe0);

        let fv2 = new_corner_vertex(k, pseudo, corner)?;
        pseudo
            .fillet_vertex_mut(fv2)?
            .set_fillet_vertex_type(FilletVertexType::RailXExtendedEdgeuse);
        let side = need(
            find_side_edgeuse(target, convex, curr)?,
            "no side edge use intersects the convex rail",
        )?;
        let rail = convex_solver.find_index_of_rail_x_side_edgeuse(target, side)?;
        pseudo.make_rail_edgeuse_on(
            convex_fg,
            rail,
            Some(fv2),
            OrientType::Unknown,
            None,
            lookup(corner, self),
            convex_solver,
            target,
        )?;
        pseudo
            .fillet_vertex_mut(fv2)?
            .set_point_class(PointClassType::Edgeuse, Some(EntityRef::Edgeuse(side)));
        self.vertices.push(fv2);

        // Make a cross section on the convex fillet
        let fe1 = pseudo.make_fillet_edge(Some(fv1), Some(fv2), Some(corner))?;
        pseudo
            .fillet_edge_mut(fe1)?
            .set_fillet_edge_type(FilletEdgeType::CrossSection);
        let primary = need(
            eu1::primary_edgeuse(pseudo.brep(), fe1)?,
            "fillet edge without primary edge use",
        )?;
        pseudo.set_edgeuse_fillet_geom(primary, Some(convex_fg))?;
        self.edges.push(fe1);

        let fv3 = new_corner_vertex(k, pseudo, corner)?;
        let record = pseudo.fillet_vertex_mut(fv3)?;
        record.set_point_class(PointClassType::Vertex, Some(EntityRef::Vertex(vertex)));
        record.set_fillet_vertex_type(FilletVertexType::OnVertex);
        self.vertices.push(fv3);

        // Find the original face
        let orig_face = face_by_two_edgeuses(target, Some(curr), Some(convex), None)?;
        let curr_edge = edge_of(target, curr)?;

        // Make an extended edge, then the second edge
        for (start, end, edge_type) in [
            (fv2, fv3, FilletEdgeType::OnExtendEdge),
            (fv3, fv0, FilletEdgeType::OnEdge),
        ] {
            let fe = pseudo.make_fillet_edge(Some(start), Some(end), Some(corner))?;
            let record = pseudo.fillet_edge_mut(fe)?;
            record.set_fillet_edge_type(edge_type);
            record.set_original_edge(Some(curr_edge));
            record.set_original_face(orig_face);
            self.edges.push(fe);
        }

        // Extend the Side Surface if necessary
        let ext_dist =
            find_surface_extension_dist_at_corner(target, convex, &mut convex_radius, vertex)?;
        let orig_face1 = need(
            face_by_two_edgeuses(target, Some(curr), Some(concave), None)?,
            "no face between the side and concave edges",
        )?;
        let surface = target
            .face(orig_face1)?
            .surface
            .as_ref()
            .ok_or_else(|| SolidError::null_pointer(OP))?;
        let new_surface = surface.create_extended_surface(k, ext_dist, ContinuityType::G1)?;
        // Copy it, when possible as an analytic surface
        let new_analytic_surface = new_surface.copy_and_add_analytics(k)?;
        FilletCorner::add_extended_original_surface_with_map(
            extended_surfaces,
            orig_face1,
            new_analytic_surface,
        );
        Ok(())
    }
}
