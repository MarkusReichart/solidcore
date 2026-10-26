//! tangency-labs/solidcore#410: a cylinder minus a screw-swept thread groove. Upstream's
//! result fails its own check: edge/vertex gaps of about 1e-5 mm against a tolerance of
//! 5e-6, sometimes "Bad Edge ValidateGeometry". Here the groove is a triangle swept
//! along a helix (`SmApiCreateCurveSweep`), cut from an r 3 shaft.
use solidcore::boolean::api_brep::boolean_difference;
use solidcore::curve::curve::Curve;
use solidcore::curve::line::Line;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::Point3;
use solidcore::primitive::api_curves::{create_helix, HelixParams};
use solidcore::primitive::api_primitives::create_cylinder;
use solidcore::primitive::primitive_creation::SweepOptions;
use solidcore::sweep::api_brep::create_curve_sweep;
use solidcore::topology::brep::p1;
use solidcore::topology::brep::p2::is_manifold_solid;
use solidcore::{Brep, Kernel};

fn volume(k: &Kernel, b: &mut Brep) -> solidcore::Result<f64> {
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    brep_volume(k, b, 1e-6, &attrs)
}

/// The largest gap between an edge's curve end and its vertex, over the larger of the
/// edge's and the vertex's tolerance, with that gap (the worst edge end).
fn worst_vertex_gap(k: &Kernel, b: &Brep) -> (f64, f64) {
    let mut worst = (0.0, 0.0);
    for e in p1::edges(b, None).unwrap() {
        let edge = b.edge(e).unwrap();
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        for (v, t) in edge
            .vertices
            .iter()
            .zip([edge.interval.min(), edge.interval.max()])
        {
            let vertex = b.vertex(*v).unwrap();
            let gap = curve.evaluate_point(k, t).unwrap().distance(vertex.point);
            let ratio = gap / edge.tolerance.max(vertex.tolerance);
            if ratio > worst.0 {
                worst = (ratio, gap);
            }
        }
    }
    worst
}

#[test]
fn issue_410_shaft_minus_thread_groove_is_valid() {
    let mut k = Kernel::default();
    let shaft = create_cylinder(&mut k, Point3::new(0.0, 0.0, 0.0), 3.0, 10.0).unwrap();
    // a pitch-1 helix on the shaft's surface from z = 1 to 9
    let helix = create_helix(&HelixParams {
        origin: Point3::new(0.0, 0.0, 1.0),
        height: 8.0,
        radius_start: 3.0,
        radius_end: 3.0,
        num_turns: 8.0,
        right_handed: true,
        tolerance: 1e-6,
    })
    .unwrap();
    let path = Curve::from(helix);
    // a triangle notch across the surface at the helix start, in the XZ plane
    let [a, b, c] = [
        Point3::new(3.3, 0.0, 0.7),
        Point3::new(2.6, 0.0, 1.0),
        Point3::new(3.3, 0.0, 1.3),
    ];
    // (each line with its NURB made: `Line::from_points` makes none, and the sweep reads
    // the profile curves' NURBs)
    let line = |p: Point3, q: Point3| {
        let mut l = Line::from_points(p, q, 3);
        l.make_nurb().unwrap();
        Curve::from(l)
    };
    let sides = [line(a, b), line(b, c), line(c, a)];
    let profile: Vec<Option<&Curve>> = sides.iter().map(Some).collect();
    let options = SweepOptions {
        cap_ends: true,
        ..SweepOptions::default()
    };
    let groove = create_curve_sweep(&mut k, &profile, Some(&path), None, None, Some(&options))
        .unwrap()
        .brep;
    // the groove on its own: a closed solid of about its Pappus volume (the triangle's
    // area 0.21 swept 8 turns around its centroid's radius 3.0667: 32.37)
    let mut groove = groove;
    let groove_manifold = is_manifold_solid(&groove).unwrap();
    let groove_volume = volume(&k, &mut groove);
    assert!(
        groove_manifold && groove_volume.as_ref().is_ok_and(|v| (v - 32.37).abs() < 1.0),
        "groove: manifold {groove_manifold}, volume {groove_volume:?}"
    );
    let threaded = boolean_difference(&mut k, shaft, groove).unwrap();
    assert!(is_manifold_solid(&threaded).unwrap(), "not a manifold solid");
    let (ratio, gap) = worst_vertex_gap(&k, &threaded);
    assert!(ratio <= 1.0, "edge/vertex gap {gap} is {ratio} times its tolerance");
}
