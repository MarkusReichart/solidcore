//! tangency-labs/solidcore#387: a box minus a cylinder given as exact rational NURBS
//! (two rational quadratic half-cylinder faces) gave SolidCraft an invalid result (edge and
//! vertex tolerances of 0.78 mm against 0.022; `ValidatePointers` fails), while the same
//! box minus SMLib's analytic cylinder is valid. Here the cylinder is made NURBS with
//! `SmBrep::TurnToNURBS`. The Rust port gives a valid cut both ways (not reproduced; see
//! `plan/findings/boolean.api_brep--issues-480-387-388-not-reproduced-regression-guards.md`);
//! this guards it.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
use solidcore::boolean::api_brep::boolean_difference;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::object::ObjectRef;
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::Point3;
use solidcore::primitive::api_primitives::{create_box, create_cylinder};
use solidcore::topology::brep::{p2::is_manifold_solid, p3};
use solidcore::topology::edge::NotifyEnv;
use solidcore::{Brep, Kernel};

fn volume(k: &Kernel, b: &mut Brep) -> solidcore::Result<f64> {
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    brep_volume(k, b, 1e-6, &attrs)
}

/// The box (-14, 0, -14)-(14, 20, 14) less the r 40 cylinder about (0, 50) along Z
/// (z -60..60), analytic or made NURBS: whether it is a manifold solid, and its volume.
fn cut(nurbs: bool) -> Result<(bool, f64), String> {
    let e = |x: solidcore::SolidError| format!("{x:?}");
    let mut k = Kernel::default();
    let block = create_box(&mut k, Point3::new(-14.0, 0.0, -14.0), 28.0, 20.0, 28.0).map_err(e)?;
    let mut tool =
        create_cylinder(&mut k, Point3::new(0.0, 50.0, -60.0), 40.0, 120.0).map_err(e)?;
    if nurbs {
        let mut records = std::collections::BTreeMap::new();
        let mut attributes = AttributeStore::default();
        let mut env = NotifyEnv {
            records: &mut records,
            attributes: &mut attributes,
            context: None,
            // caller keys below the kernel's reserved range (as patterned_cut.rs names them)
            key_of: &|x| {
                let r = p3::topology_ref(0, x);
                ObjectRef {
                    type_code: r.type_code,
                    key: ObjectKey(r.type_code as u64 * 1_000_000 + (r.key.0 & 0xffff_ffff)),
                }
            },
        };
        p3::turn_to_nurbs(&mut k, &mut tool, &mut env).map_err(e)?;
    }
    let mut out = boolean_difference(&mut k, block, tool).map_err(e)?;
    solidcore::topology::integrity::check(&out).map_err(e)?;
    let solid = is_manifold_solid(&out).map_err(e)?;
    Ok((solid, volume(&k, &mut out).map_err(e)?))
}

/// The box less its part inside the cylinder: 28 (28 * 20 - A), A the area of the part
/// of [-14, 14] x [0, 20] within 40 of (0, 50).
fn want() -> f64 {
    let s = (1600.0_f64 - 196.0).sqrt();
    let area = 2.0 * (7.0 * s + 800.0 * 0.35_f64.asin()) - 30.0 * 28.0;
    28.0 * (28.0 * 20.0 - area)
}

#[test]
fn issue_387_box_minus_nurbs_cylinder_matches_the_analytic_cut() {
    for nurbs in [false, true] {
        let (solid, v) = cut(nurbs).unwrap_or_else(|m| panic!("nurbs {nurbs}: {m}"));
        assert!(solid, "nurbs {nurbs}: not a manifold solid");
        assert!(
            (v - want()).abs() <= 1e-5 * want(),
            "nurbs {nurbs}: {v} vs {}",
            want()
        );
    }
}
