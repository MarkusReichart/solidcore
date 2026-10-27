//! Probes for SolidCraft fillet cases that stop at `MergeCurveOnSurfaces`
//! (`SM_ERR_NULL_POINTER`): a single box-edge round, varied by box size and edge, to
//! tell whether the failure follows the plate's size or the edge chosen. The box-edge
//! unit test (`fillet_one_box_edge_end_to_end`: a 10 x 20 x 30 box, its first edge)
//! passes; SolidCraft #2707's `round_on_plate_100` (a 100 x 100 x 10 plate, the top
//! edge at y = 100) does not.

use solidcore::blend::fillet_executive::api::fillet_edges;
use solidcore::kernel::Kernel;
use solidcore::math::Point3;
use solidcore::primitive::api_primitives::create_box;
use solidcore::topology::{edge, Brep, EdgeId};

/// The edge of `b` whose two ends satisfy `on`.
fn edge_where(k: &Kernel, b: &Brep, on: impl Fn(Point3) -> bool) -> EdgeId {
    let found: Vec<EdgeId> = b
        .edges()
        .ids()
        .filter(|&e| {
            let curve = edge::curve(b, e).unwrap().unwrap();
            let ivl = edge::interval(b, e).unwrap();
            on(curve.evaluate_point(k, ivl.min()).unwrap())
                && on(curve.evaluate_point(k, ivl.max()).unwrap())
        })
        .collect();
    assert_eq!(found.len(), 1, "edges matching: {found:?}");
    found[0]
}

fn round(size: [f64; 3], pick: impl Fn(Point3) -> bool) {
    let mut k = Kernel::default();
    let mut b = create_box(
        &mut k,
        Point3::new(0.0, 0.0, 0.0),
        size[0],
        size[1],
        size[2],
    )
    .unwrap();
    let e = edge_where(&k, &b, pick);
    fillet_edges(&mut k, &mut b, &[e], 1.0, 1, 1, 1.0).unwrap();
    let faces = solidcore::topology::brep::p1::faces(&b, None).unwrap().len();
    assert_eq!(faces, 7);
}

/// SolidCraft #2707 `round_on_plate_100`: the plate's top edge at y = 100.
#[test]
fn plate_100_top_edge_at_y_100() {
    round([100.0, 100.0, 10.0], |p| p.y == 100.0 && p.z == 10.0);
}

/// The same plate, the edge the unit test takes (the first edge of the box).
#[test]
fn plate_100_first_edge() {
    let mut k = Kernel::default();
    let mut b = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 100.0, 100.0, 10.0).unwrap();
    let e = b.edges().ids().next().unwrap();
    fillet_edges(&mut k, &mut b, &[e], 1.0, 1, 1, 1.0).unwrap();
    let faces = solidcore::topology::brep::p1::faces(&b, None).unwrap().len();
    assert_eq!(faces, 7);
}

/// The unit test's box, the edge in the same place as #2707's (top, far side).
#[test]
fn box_10_20_30_top_edge_at_far_side() {
    round([10.0, 20.0, 30.0], |p| p.y == 20.0 && p.z == 30.0);
}

/// The unit test's box and edge, through this entry, as the control.
#[test]
fn box_10_20_30_first_edge_control() {
    let mut k = Kernel::default();
    let mut b = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 10.0, 20.0, 30.0).unwrap();
    let e = b.edges().ids().next().unwrap();
    fillet_edges(&mut k, &mut b, &[e], 1.0, 1, 1, 1.0).unwrap();
    let faces = solidcore::topology::brep::p1::faces(&b, None).unwrap().len();
    assert_eq!(faces, 7);
}
