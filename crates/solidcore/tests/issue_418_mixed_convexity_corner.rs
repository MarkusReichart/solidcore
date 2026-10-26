//! The mixed-convexity Nx1 corner (`SmFilletNx1Corner::MakeCornerTopology`'s capping
//! branch, `sm_FindSurfaceExtensionDistAtCorner`; tangency-labs/solidcore#418): a round
//! on a convex edge of an L-shaped block that ends where the other side edge at that
//! vertex is the concave inner edge. The corner builds its end rib, the end face's
//! extension and, at the corner's vertex, a fillet vertex and the two fillet edges to the
//! rails' ends, instead of stopping `Unported`.

use solidcore::blend::fillet_executive::session::{FilletObject, FilletSession};
use solidcore::boolean::api_brep::boolean_union;
use solidcore::math::Point3;
use solidcore::primitive::api_primitives::create_box;
use solidcore::topology::{Brep, EdgeId, VertexId};
use solidcore::Kernel;

/// The vertex of `b` at `p`.
fn vertex_at(b: &Brep, p: Point3) -> VertexId {
    b.vertices()
        .iter()
        .find(|(_, v)| v.point.distance(p) < 1e-9)
        .map(|(id, _)| id)
        .unwrap()
}

/// The edge of `b` between the vertices at `p` and `q`.
fn edge_between(b: &Brep, p: Point3, q: Point3) -> EdgeId {
    let (a, z) = (vertex_at(b, p), vertex_at(b, q));
    b.edges()
        .iter()
        .find(|(_, e)| {
            (e.vertices[0] == a && e.vertices[1] == z) || (e.vertices[0] == z && e.vertices[1] == a)
        })
        .map(|(id, _)| id)
        .unwrap()
}

#[test]
fn issue_418_mixed_convexity_nx1_corner_builds_its_cap_topology() {
    let mut k = Kernel::default();
    // an L-shaped block: [0,20]x[0,10]x[0,10] united with [0,10]x[0,20]x[0,10]; its
    // inner vertical edge at x = y = 10 is concave
    let a = create_box(&mut k, Point3::ZERO, 20.0, 10.0, 10.0).unwrap();
    let b = create_box(&mut k, Point3::ZERO, 10.0, 20.0, 10.0).unwrap();
    let mut target = boolean_union(&mut k, a, b).unwrap();
    let corner_point = Point3::new(10.0, 10.0, 10.0);
    let edge = edge_between(&target, Point3::new(20.0, 10.0, 10.0), corner_point);
    let corner_vertex = vertex_at(&target, corner_point);

    let mut session = FilletSession::new(&mut k, &mut target).unwrap();
    session
        .prepare_constant_radius_edges(&mut k, &[edge], 1.0, 1, 1, 1.0)
        .unwrap();
    session.create_corners(&k, None).unwrap();
    let id = session
        .executive()
        .corners
        .iter()
        .find_map(|corner| match corner {
            FilletObject::Corner(id) => {
                (session.corner(*id).unwrap().vertex == Some(corner_vertex)).then_some(*id)
            }
            _ => None,
        })
        .expect("a corner at the concave edge's top vertex");
    session.make_nx1_corner_topology(&k, id).unwrap();
    let corner = session.corner(id).unwrap();
    // the two rail vertices and the fillet vertex on the corner's own vertex; the end
    // rib and the two fillet edges from the rails' ends to that vertex
    assert_eq!(corner.vertices.len(), 3);
    assert_eq!(corner.edges.len(), 3);
}
