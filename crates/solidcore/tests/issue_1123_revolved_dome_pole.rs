//! tangency-labs/solidcore#1123 (SolidCraft rig shell/4847, step 6): a quarter disc of
//! radius 20 revolved a whole turn, its arc a sketch's arc through (20, 0), a rounded
//! (14.142136, 14.142136) and (0, 20), so centred 1.28e-6 off the axis with a radius
//! 1.28e-6 short of 20. In the SolidCraft fork the dome came out with a pole row of
//! radius 1.26e-6 and an off-axis pole vertex: a hole. That came from the fork's analytic
//! sphere branch (`SmRotationalSweepGeometry::AnalyticSurface`, not in `c9c5979`).
//! Upstream, and the port, sweep each curve of a whole turn into a surface of revolution
//! (`SmPrimitiveCreation::CreateRotationalSweep`), whose row at the arc's end on the axis
//! is that end itself. This checks the port's result has no hole there.
//!
//! The profile is drawn in the XY plane and revolved about Y (the issue's XZ profile
//! about Z, turned): (0,0)→(20,0), the arc, (0,20)→(0,0).

#![allow(clippy::unwrap_used, clippy::indexing_slicing, clippy::print_stdout)]

use solidcore::curve::curve::Curve;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::kernel::Kernel;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::{Point3, Vec3};
use solidcore::primitive::api_curves::{create_arc, create_line_segment};
use solidcore::sweep::api_brep::create_rotational_sweep;
use solidcore::topology::brep::p1 as bl;
use solidcore::topology::brep::p2::is_manifold_solid;

/// The arc's centre offset from both axes and its radius (the issue's numbers).
const C: f64 = 1.2846627950203038e-06;
const R: f64 = 19.99999871533725;
/// Its start angle (radians) and sweep: from (20, 0) to (0, 20).
const A0: f64 = -6.42331438769116e-08;
const SWEEP: f64 = 1.5707964552611844;

#[test]
fn issue_1123_revolved_dome_has_its_pole_on_the_axis() {
    let mut k = Kernel::default();
    let to_deg = 180.0 / std::f64::consts::PI;
    let bottom = Curve::from(
        create_line_segment(Point3::new(0.0, 0.0, 0.0), Point3::new(20.0, 0.0, 0.0)).unwrap(),
    );
    let arc = Curve::from(
        create_arc(Point3::new(C, C, 0.0), R, A0 * to_deg, (A0 + SWEEP) * to_deg).unwrap(),
    );
    let side = Curve::from(
        create_line_segment(Point3::new(0.0, 20.0, 0.0), Point3::new(0.0, 0.0, 0.0)).unwrap(),
    );
    let mut b = create_rotational_sweep(
        &mut k,
        &[Some(&bottom), Some(&arc), Some(&side)],
        Point3::ZERO,
        Vec3::new(0.0, 1.0, 0.0),
        360.0,
        true,
    )
    .unwrap();
    // the integrity check, reported now and asserted last, so the pole and volume
    // checks still run; on failure print each face's `shell` field, the shell of its
    // face uses and the live shells with their face lists
    let integrity = solidcore::topology::integrity::check(&b);
    if integrity.is_err() {
        let live: Vec<_> = b.shells().iter().map(|(id, s)| (id, s.faces.clone())).collect();
        println!("live shells (id, faces): {live:?}");
        for (f, face) in b.faces().iter() {
            let uses = face.faceuse.map(|fu| {
                let up = b.faceuse(fu).unwrap();
                (fu, up.shell, up.mate.map(|m| b.faceuse(m).unwrap().shell))
            });
            println!(
                "face {f:?}: shell field {:?}, faceuse (id, shell, mate shell) {uses:?}",
                face.shell
            );
        }
    }
    assert!(is_manifold_solid(&b).unwrap(), "not a manifold solid");

    // the pole: every vertex near the top lies on the axis
    for v in bl::vertices(&b, None).unwrap() {
        let p = b.vertex(v).unwrap().point;
        if p.y > 19.0 {
            let off_axis = (p.x * p.x + p.z * p.z).sqrt();
            assert!(off_axis < 1e-9, "pole vertex {p:?} is {off_axis:e} off the axis");
        }
    }

    // a hemisphere of radius 20 to the profile's 1.3e-6 (a hole of the issue's size,
    // a pole parallel of radius 1.26e-6, would not show in the volume; the vertex and
    // the solid check above are the real test)
    let store = AttributeStore::new();
    let key = |_| ObjectKey(0);
    let attributes = PropertyAttributes {
        store: &store,
        key_of: &key,
    };
    let volume = brep_volume(&k, &mut b, 1e-9, &attributes).unwrap();
    let expected = 2.0 / 3.0 * std::f64::consts::PI * 20.0_f64.powi(3);
    assert!(
        ((volume - expected) / expected).abs() < 1e-5,
        "volume {volume} vs {expected}"
    );
    integrity.unwrap();
}
