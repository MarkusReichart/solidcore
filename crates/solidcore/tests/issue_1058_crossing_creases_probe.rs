//! tangency-labs/solidcore#1058 (WR-306) probe, report only: two convex creases crossing
//! with a tiny overlap. In the SolidCraft fork's rig case, creases passing about 1e-5
//! apart gave an intersection curve 1.4e-5 long with both ends in one proto vertex on
//! one vertex definition; upstream kept it as a closed edge on that vertex (it is longer
//! than twice its own tolerance) and imprinted a speck face. With
//! `Kernel::fixes.degenerate_proto_edge_in_one_vertex` it is that vertex.
//!
//! Here: the box A = [0,10]^3 and a box B whose edge (its crease) crosses A's edge
//! y = z = 10 at right angles, outside A except for a depth `d` along A's edge's outward
//! diagonal. The overlap is a tetrahedron of size about `d`, so the union is 1000 +
//! B's volume less about d^3. Each depth runs with the flag off and on; the PROBE lines
//! say which reach the branch (a face or edge count, or a check, that differs).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers

use solidcore::boolean::api_brep::boolean_union;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::{Point3, Vec3};
use solidcore::primitive::api_primitives::create_box;
use solidcore::topology::brep::p1;
use solidcore::topology::brep::p2::is_manifold_solid;
use solidcore::topology::proto_topology::p1::api_transform::{
    api_rotate, api_translate, TransformObject,
};
use solidcore::{Brep, Kernel};

fn volume(k: &Kernel, b: &mut Brep) -> solidcore::Result<f64> {
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    brep_volume(k, b, 1e-8, &attrs)
}

/// The largest gap between an edge's curve end and its vertex, and how many edge ends
/// have no live vertex. `topological` reads the ends through the edge uses
/// (`GetStartVertex`/`GetEndVertex`, `edge::start_end_vertices`); otherwise from the
/// `Edge::vertices` field, kept for comparison (the field is a co-ported shortcut that
/// the topology edits may not maintain).
fn worst_vertex_gap(k: &Kernel, b: &Brep, topological: bool) -> (f64, usize) {
    let mut worst = 0.0_f64;
    let mut dangling = 0;
    for e in p1::edges(b, None).unwrap() {
        let edge = b.edge(e).unwrap();
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let ends = if topological {
            let (s, t) = solidcore::topology::edge::start_end_vertices(b, e).unwrap();
            [s, t]
        } else {
            [Some(edge.vertices[0]), Some(edge.vertices[1])]
        };
        for (v, t) in ends
            .into_iter()
            .zip([edge.interval.min(), edge.interval.max()])
        {
            let Some(Ok(vertex)) = v.map(|v| b.vertex(v)) else {
                dangling += 1;
                continue;
            };
            worst = worst.max(curve.evaluate_point(k, t).unwrap().distance(vertex.point));
        }
    }
    (worst, dangling)
}

/// B: a 10 x 10 x 20 box whose edge (x = y = 0, z -10..10) is turned onto the line
/// through Q along (0, -1, 1), its body on the side (0, 1, 1) of it.
fn crease(k: &mut Kernel, d: f64) -> Brep {
    let mut b = create_box(k, Point3::new(0.0, 0.0, -10.0), 10.0, 10.0, 20.0).unwrap();
    let o = Point3::new(0.0, 0.0, 0.0);
    // its bisector (1, 1, 0) onto +y, then +y onto (0, 1, 1) and +z onto (0, -1, 1)
    api_rotate(
        k,
        Some(TransformObject::Brep(&mut b)),
        o,
        Vec3::new(0.0, 0.0, 1.0),
        45.0,
    )
    .unwrap();
    api_rotate(
        k,
        Some(TransformObject::Brep(&mut b)),
        o,
        Vec3::new(1.0, 0.0, 0.0),
        45.0,
    )
    .unwrap();
    let s = d / 2.0_f64.sqrt();
    api_translate(
        k,
        Some(TransformObject::Brep(&mut b)),
        Vec3::new(5.0, 10.0 - s, 10.0 - s),
    )
    .unwrap();
    b
}

fn run(d: f64, fix: bool) {
    let mut k = Kernel::default();
    k.fixes.degenerate_proto_edge_in_one_vertex = fix;
    let a = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 10.0, 10.0, 10.0).unwrap();
    let b = crease(&mut k, d);
    match boolean_union(&mut k, a, b) {
        Err(e) => eprintln!("PROBE issue_1058 d={d:e} fix={fix}: Err({e:?})"),
        Ok(mut u) => {
            let faces = p1::faces(&u, None).unwrap().len();
            let edges = p1::edges(&u, None).unwrap().len();
            let vertices = p1::vertices(&u, None).unwrap().len();
            let manifold = is_manifold_solid(&u);
            let integrity = solidcore::topology::integrity::check(&u);
            let (gap, dangling) = worst_vertex_gap(&k, &u, true);
            let (field_gap, field_dangling) = worst_vertex_gap(&k, &u, false);
            let v = volume(&k, &mut u);
            eprintln!(
                "PROBE issue_1058 d={d:e} fix={fix}: Ok; faces={faces} edges={edges} vertices={vertices} manifold={manifold:?} worst_vertex_gap={gap:e} dangling_vertex_refs={dangling} integrity={integrity:?} (Edge::vertices field: gap={field_gap:e} dangling={field_dangling}) volume={v:?} (~3000)"
            );
        }
    }
}

#[test]
fn issue_1058_crossing_creases_probe() {
    for d in [3e-6, 1e-5, 3e-5, 1e-4] {
        for fix in [false, true] {
            run(d, fix);
        }
    }
}

/// The measure on a trivial union of two overlapping axis-aligned boxes: through the
/// edge uses it must find every end on a live vertex with a gap at round-off.
#[test]
fn issue_1058_measure_control_two_boxes() {
    let mut k = Kernel::default();
    let a = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 10.0, 10.0, 10.0).unwrap();
    let b = create_box(&mut k, Point3::new(5.0, 5.0, 5.0), 10.0, 10.0, 10.0).unwrap();
    let u = boolean_union(&mut k, a, b).unwrap();
    solidcore::topology::integrity::check(&u).unwrap();
    let (gap, dangling) = worst_vertex_gap(&k, &u, true);
    let (field_gap, field_dangling) = worst_vertex_gap(&k, &u, false);
    eprintln!(
        "PROBE measure control: edge uses gap={gap:e} dangling={dangling}; Edge::vertices field gap={field_gap:e} dangling={field_dangling}"
    );
    assert_eq!(dangling, 0);
    assert!(gap < 1e-9, "gap {gap:e}");
}
