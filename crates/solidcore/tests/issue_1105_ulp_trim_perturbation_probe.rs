//! tangency-labs/solidcore#1105 probe (hypothesis check, Linux): on Windows GNU the pole-free
//! sphere of radius 5 at (1, -2, 0.5), construction tolerance 1e-9, misses its analytic
//! area by 5.5e-9 relative at request accuracy 1e-10. Nearly all of it is in `FaceId(3)`,
//! and the only reported difference between the hosts is about one ulp in a few fitted UV
//! trim control coordinates. The hypothesis is that the whole-trim adaptive integration
//! at a request this tight is sensitive enough that a one-ulp change in a trim is enough.
//!
//! This builds the same sphere on this host, nudges every UV trim control coordinate of one
//! face by one ulp (up or down), and reruns `brep_properties` at 1e-10. If an area error
//! near 5e-9 appears, the sensitivity is the integrator's and not Windows'. It reports only
//! (run with `--ignored --nocapture`); nothing is asserted.
//!
//! The nudge happens after construction, so caches the construction may have built for
//! the face are not rebuilt; the printed per-face areas show whether the nudge reached the
//! integration (a face whose area is bit-identical to the baseline was not affected).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
use solidcore::curve::curve::CurveKind;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::mass::integrator::properties::{
    brep_properties, MassPropertyBehavior, PropertyAttributes,
};
use solidcore::math::{Point3, Vec3};
use solidcore::nurbs::b_spline_curve::ControlPointForm;
use solidcore::primitive::api_primitives::create_sphere_no_pole_with_tolerance;
use solidcore::query::brep::{calculate_bounding_box, BoundingBoxOptions};
use solidcore::topology::edgeuse::p1 as eu;
use solidcore::topology::FaceId;
use solidcore::{Brep, Kernel};

const RADIUS: f64 = 5.0;
const ACCURACY: f64 = 1e-10;

fn sphere(k: &mut Kernel) -> Brep {
    create_sphere_no_pole_with_tolerance(k, Vec3::new(1.0, -2.0, 0.5), RADIUS, 1e-9).unwrap()
}

fn area(k: &Kernel, b: &mut Brep) -> f64 {
    let bbox = calculate_bounding_box(k, &*b, BoundingBoxOptions::default()).unwrap();
    let store = AttributeStore::default();
    let key_of = |_| ObjectKey(0);
    brep_properties(
        k,
        b,
        ACCURACY,
        bbox.mid(),
        MassPropertyBehavior::default(),
        &PropertyAttributes {
            store: &store,
            key_of: &key_of,
        },
    )
    .unwrap()
    .area
}

/// Moves every UV trim control coordinate (x and y, stored form) of `face`'s edge uses one
/// ulp up (`up`) or down; returns how many coordinates moved.
fn nudge_face_trims(b: &mut Brep, face: FaceId, up: bool) -> usize {
    let edgeuses: Vec<_> = b
        .edgeuses()
        .ids()
        .filter(|&e| eu::face(b, e).ok().flatten() == Some(face))
        .collect();
    let mut moved = 0;
    for e in edgeuses {
        let Some(curve) = b.edgeuse_mut(e).unwrap().uv_trim_curve.as_mut() else {
            continue;
        };
        let CurveKind::BSpline(c) = curve.kind_mut() else {
            continue;
        };
        let points = c.control_points().unwrap().to_vec();
        for (i, p) in points.iter().enumerate() {
            let step = |v: f64| if up { v.next_up() } else { v.next_down() };
            c.set_control_point(
                ControlPointForm::HomogeneousRational,
                i,
                Point3::new(step(p.x), step(p.y), p.z),
                p.w,
            )
            .unwrap();
            moved += 2;
        }
    }
    moved
}

#[test]
#[ignore = "probe for #1105: reports only"]
fn issue_1105_one_ulp_trim_nudge_probe() {
    let want = 4.0 * std::f64::consts::PI * RADIUS * RADIUS;
    let mut k = Kernel::default();
    let mut base = sphere(&mut k);
    let a0 = area(&k, &mut base);
    println!(
        "baseline area {a0:.15} rel {:.3e}",
        (a0 - want) / want
    );
    let faces: Vec<FaceId> = base.faces().ids().collect();
    for face in faces {
        for up in [true, false] {
            let mut k = Kernel::default();
            let mut b = sphere(&mut k);
            let moved = nudge_face_trims(&mut b, face, up);
            let a = area(&k, &mut b);
            println!(
                "face {face:?} {} ({moved} coords): area {a:.15} rel {:.3e} (vs baseline {:.3e})",
                if up { "+1 ulp" } else { "-1 ulp" },
                (a - want) / want,
                (a - a0) / want
            );
        }
    }
}
