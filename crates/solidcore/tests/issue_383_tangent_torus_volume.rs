//! tangency-labs/solidcore#383, the torus rows (seeds 4 and 86): a box and a torus whose
//! tube is tangent to one of the box's faces along a circle arc inside that face. The
//! intersection and the cut are checked against a reference volume of box n torus from a
//! deterministic midpoint grid (the torus's inside test), so either can be judged on its
//! own.
use solidcore::boolean::api_brep::{boolean_difference, boolean_intersection};
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::Point3;
use solidcore::primitive::api_primitives::{create_box, create_torus};
use solidcore::{Brep, Kernel};

fn volume(k: &Kernel, b: &mut Brep) -> solidcore::Result<f64> {
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    brep_volume(k, b, 1e-6, &attrs)
}

struct Case {
    name: &'static str,
    lo: [f64; 3],
    hi: [f64; 3],
    centre: [f64; 3],
    big: f64,
    small: f64,
}

/// The volume of box n torus (torus axis +Z) by the midpoint rule on an `n`^3 grid over
/// the box.
fn reference(c: &Case, n: usize) -> f64 {
    let step = [0, 1, 2].map(|i| (c.hi[i] - c.lo[i]) / n as f64);
    let mut inside = 0usize;
    for i in 0..n {
        let x = c.lo[0] + (i as f64 + 0.5) * step[0] - c.centre[0];
        for j in 0..n {
            let y = c.lo[1] + (j as f64 + 0.5) * step[1] - c.centre[1];
            let ring = (x * x + y * y).sqrt() - c.big;
            for l in 0..n {
                let z = c.lo[2] + (l as f64 + 0.5) * step[2] - c.centre[2];
                if ring * ring + z * z < c.small * c.small {
                    inside += 1;
                }
            }
        }
    }
    inside as f64 * step[0] * step[1] * step[2]
}

fn operands(k: &mut Kernel, c: &Case) -> (Brep, Brep) {
    let size = [0, 1, 2].map(|i| c.hi[i] - c.lo[i]);
    let block = create_box(k, Point3::new(c.lo[0], c.lo[1], c.lo[2]), size[0], size[1], size[2])
        .unwrap();
    let ring = create_torus(
        k,
        Point3::new(c.centre[0], c.centre[1], c.centre[2]),
        c.big,
        c.small,
    )
    .unwrap();
    (block, ring)
}

#[test]
fn issue_383_box_and_tangent_torus_volumes() {
    let cases = [
        Case {
            name: "seed 4 (tube tangent to the box's bottom z = 22.5)",
            lo: [22.5, 25.0, 22.5],
            hi: [30.0, 32.5, 50.0],
            centre: [20.0, 22.5, 25.0],
            big: 10.0,
            small: 2.5,
        },
        Case {
            name: "seed 86 (tube tangent to the box's top z = 32.5)",
            lo: [30.0, 20.0, 15.0],
            hi: [35.0, 30.0, 32.5],
            centre: [20.0, 17.5, 30.0],
            big: 12.5,
            small: 2.5,
        },
    ];
    let mut report = Vec::new();
    let mut failed = 0;
    for c in &cases {
        let want = reference(c, 240);
        let block_volume: f64 = (0..3).map(|i| c.hi[i] - c.lo[i]).product();
        let tol = 2e-3 * block_volume;

        let mut k = Kernel::default();
        let (a, b) = operands(&mut k, c);
        let inter = boolean_intersection(&mut k, a, b).and_then(|mut r| volume(&k, &mut r));
        let (a, b) = operands(&mut k, c);
        let cut = boolean_difference(&mut k, a, b).and_then(|mut r| volume(&k, &mut r));

        let inter_ok = inter.as_ref().is_ok_and(|v| (v - want).abs() <= tol);
        let cut_ok = cut
            .as_ref()
            .is_ok_and(|v| (v - (block_volume - want)).abs() <= tol);
        failed += usize::from(!inter_ok) + usize::from(!cut_ok);
        report.push(format!(
            "{}: reference box n torus {want:.3}; intersection {inter:?} ({}); \
             cut {cut:?}, want {:.3} ({})",
            c.name,
            if inter_ok { "ok" } else { "WRONG" },
            block_volume - want,
            if cut_ok { "ok" } else { "WRONG" },
        ));
    }
    assert_eq!(failed, 0, "{report:#?}");
}
