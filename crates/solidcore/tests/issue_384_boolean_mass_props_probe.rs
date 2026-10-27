//! tangency-labs/solidcore#384: SolidCraft saw SMLib's exact volume come out wrong on some
//! box/cylinder Boolean results whose bodies were right (a cut larger than the body it was
//! cut from). The fork did not reproduce it with SMLib's own cylinders, and neither does the
//! port. This regression test builds the issue's three described seeds (2, 6, 15) with the
//! port's primitives, runs the union, intersection and cut, and measures each result with
//! `brep_volume` against independent checks:
//! - the closed form of seed 2's intersection, the issue's 2688.28;
//! - inclusion-exclusion, V(A ∪ B) + V(A ∩ B) = V(A) + V(B), with the operands' closed-form
//!   volumes;
//! - V(A − B) = V(A) − V(A ∩ B).
//!
//! At the integrator accuracy 1e-6 every identity holds within 2e-6 relative, and no cut
//! exceeds its operand. (At 1e-8, beyond `SmApiBrepComputeMassProperties`' documented
//! [1e-4, 1e-1], the residuals reached 1.4e-7: the error control does not meet a request
//! that tight, so it is not asserted.)
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
use core::f64::consts::PI;
use solidcore::boolean::api_brep::{boolean_difference, boolean_intersection, boolean_union};
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::{Point3, Vec3};
use solidcore::primitive::api_primitives::{create_box, create_cylinder};
use solidcore::topology::proto_topology::p1::api_transform::{
    api_rotate, api_translate, TransformObject,
};
use solidcore::{Brep, Kernel};

/// The accuracy asked of the integrator.
const ACCURACY: f64 = 1e-6;
/// The relative bound on each identity's residual.
const BOUND: f64 = 2e-6;

fn volume(k: &Kernel, b: &mut Brep, accuracy: f64) -> solidcore::Result<f64> {
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    brep_volume(k, b, accuracy, &attrs)
}

#[derive(Clone, Copy, Debug)]
enum Solid {
    /// Corners `a` and `b`.
    Box([f64; 3], [f64; 3]),
    /// Axis (0 X, 1 Y, 2 Z), base, radius, height.
    Cyl(usize, [f64; 3], f64, f64),
}

impl Solid {
    fn exact_volume(self) -> f64 {
        match self {
            Solid::Box(a, b) => (b[0] - a[0]) * (b[1] - a[1]) * (b[2] - a[2]),
            Solid::Cyl(_, _, r, h) => PI * r * r * h,
        }
    }

    fn build(self, k: &mut Kernel) -> solidcore::Result<Brep> {
        match self {
            Solid::Box(a, b) => create_box(
                k,
                Point3::new(a[0], a[1], a[2]),
                b[0] - a[0],
                b[1] - a[1],
                b[2] - a[2],
            ),
            Solid::Cyl(axis, base, r, h) => {
                let mut c = create_cylinder(k, Point3::new(0.0, 0.0, 0.0), r, h)?;
                let o = Point3::new(0.0, 0.0, 0.0);
                match axis {
                    // +Z turned onto +X about Y, onto +Y about X
                    0 => api_rotate(k, Some(TransformObject::Brep(&mut c)), o, Vec3::Y, 90.0)?,
                    1 => api_rotate(k, Some(TransformObject::Brep(&mut c)), o, Vec3::X, -90.0)?,
                    _ => {}
                }
                api_translate(
                    k,
                    Some(TransformObject::Brep(&mut c)),
                    Vec3::new(base[0], base[1], base[2]),
                )?;
                Ok(c)
            }
        }
    }
}

#[derive(Clone, Copy, Debug)]
enum Op {
    Union,
    Intersect,
    Cut,
}

fn run(a: Solid, b: Solid, op: Op, accuracy: f64) -> solidcore::Result<f64> {
    let mut k = Kernel::default();
    let x = a.build(&mut k)?;
    let y = b.build(&mut k)?;
    let mut r = match op {
        Op::Union => boolean_union(&mut k, x, y)?,
        Op::Intersect => boolean_intersection(&mut k, x, y)?,
        Op::Cut => boolean_difference(&mut k, x, y)?,
    };
    volume(&k, &mut r, accuracy)
}

/// The issue's seed 2 intersection: the r 10 Y cylinder's disc about (27.5, 15) in xz
/// clipped to x 20..35 (two segments 7.5 from the centre cut off; z 5..25 lies inside
/// 2.5..25), over the y overlap 12.5..22.5.
fn seed2_intersection() -> f64 {
    let (r, d): (f64, f64) = (10.0, 7.5);
    let segment = r * r * (d / r).acos() - d * (r * r - d * d).sqrt();
    (PI * r * r - 2.0 * segment) * 10.0
}

#[test]
fn issue_384_boolean_volumes_hold_their_identities() {
    let seeds: [(u32, Solid, Solid); 3] = [
        (
            2,
            Solid::Cyl(1, [27.5, 10.0, 15.0], 10.0, 12.5),
            Solid::Box([20.0, 12.5, 2.5], [35.0, 37.5, 25.0]),
        ),
        (
            6,
            Solid::Cyl(1, [10.0, 7.5, 12.5], 10.0, 27.5),
            Solid::Cyl(2, [17.5, 5.0, 20.0], 5.0, 15.0),
        ),
        (
            15,
            Solid::Cyl(0, [0.0, 20.0, 5.0], 7.5, 25.0),
            Solid::Box([0.0, 22.5, 7.5], [15.0, 32.5, 20.0]),
        ),
    ];
    for (seed, a, b) in seeds {
        let (va, vb) = (a.exact_volume(), b.exact_volume());
        let u = run(a, b, Op::Union, ACCURACY).unwrap();
        let i = run(a, b, Op::Intersect, ACCURACY).unwrap();
        let c = run(a, b, Op::Cut, ACCURACY).unwrap();
        let union_residual = (u + i - va - vb) / (va + vb);
        assert!(
            union_residual.abs() <= BOUND,
            "seed {seed}: (U + I - A - B) / (A + B) = {union_residual:e} (U {u}, I {i})"
        );
        let cut_residual = (c - (va - i)) / va;
        assert!(
            cut_residual.abs() <= BOUND,
            "seed {seed}: (C - (A - I)) / A = {cut_residual:e} (C {c}, I {i})"
        );
        assert!(c <= va, "seed {seed}: cut {c} exceeds its operand {va}");
        if seed == 2 {
            let want = seed2_intersection();
            assert!(
                ((i - want) / want).abs() <= BOUND,
                "seed 2: intersection {i}, closed form {want}"
            );
        }
    }
}
