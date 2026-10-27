//! tangency-labs/solidcore#388: a 20 cube joined with a hollow cube (a 50 cube less a 40
//! cube, the 20 cube inside the cavity), then cut by a 10 cube. SolidCraft reports that the
//! result fails its own check (`SmBrep` test #12, "Inner shell is not geometrically inside
//! its containing `outerShell`"). The issue does not place the 10 cube; here it overlaps the
//! island's corner (10..20 against the island's 15..35). The Rust port's result is valid
//! (not reproduced; see
//! `plan/findings/boolean.api_brep--issues-480-387-388-not-reproduced-regression-guards.md`);
//! this guards it.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
use solidcore::boolean::api_brep::{boolean_difference, boolean_union};
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::Point3;
use solidcore::primitive::api_primitives::create_box;
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

fn cube(k: &mut Kernel, at: f64, size: f64) -> Brep {
    create_box(k, Point3::new(at, at, at), size, size, size).unwrap()
}

#[test]
fn issue_388_island_in_a_void_stays_valid() {
    let mut k = Kernel::default();
    let (outer, inner) = (cube(&mut k, 0.0, 50.0), cube(&mut k, 5.0, 40.0));
    let shell = boolean_difference(&mut k, outer, inner).unwrap();
    let island = cube(&mut k, 15.0, 20.0);
    let joined = boolean_union(&mut k, island, shell).unwrap();
    let tool = cube(&mut k, 10.0, 10.0);
    let mut cut = boolean_difference(&mut k, joined, tool).unwrap();
    solidcore::topology::integrity::check(&cut).unwrap();
    assert!(is_manifold_solid(&cut).unwrap());
    let want = 50.0_f64.powi(3) - 40.0_f64.powi(3) + 20.0_f64.powi(3) - 5.0_f64.powi(3);
    let got = volume(&k, &mut cut).unwrap();
    assert!((got - want).abs() <= 1e-6 * want, "{got} vs {want}");
}
