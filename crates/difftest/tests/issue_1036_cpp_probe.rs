//! tangency-labs/solidcore#1036 on the C++ reference (needs `--features cpp`): the two
//! capped degree-1 skins of `crates/solidcore/tests/issue_1036_coincident_skin_sides.rs`
//! (A through sections 0..4, B through 2..6), built by `SmApiCreateSkinPrimitive` from the
//! exact NURB section curves the Rust test builds, then A − B, A ∩ B and A ∪ B through the
//! SM_API Booleans. It prints the C++ and the Rust results side by side; it asserts only
//! that the C++ side ran. One-off diagnostic for the #1036 record.
#![cfg(feature = "cpp")]
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::print_stdout
)]

use solidcore::boolean::api_brep::{boolean_difference, boolean_intersection, boolean_union};
use solidcore::curve::curve::Curve;
use solidcore::math::Point3;
use solidcore::nurbs::b_spline_curve::BSplineCurve;
use solidcore::nurbs::crv_fit::p1::ParamMethod;
use solidcore::sweep::api_surfaces::create_skin_primitive;
use solidcore::topology::Brep;
use solidcore::Kernel;
use solidcore_ffi::units::probe_issue_1036::{self as ffi, Nurb, Op, Skin};

/// Section `i` of the #1036 test: a closed cubic interpolant through 12 points of an
/// ellipse at height `10 i`.
fn section(i: usize) -> BSplineCurve {
    let i = f64::from(u32::try_from(i).unwrap());
    let (a, b, z) = (20.0 * (1.0 - 0.04 * i), 12.0 * (1.0 - 0.03 * i), 10.0 * i);
    let mut pts: Vec<Point3> = (0..12)
        .map(|j| {
            let t = 2.0 * std::f64::consts::PI * f64::from(j) / 12.0;
            Point3::new(a * t.cos(), b * t.sin(), z)
        })
        .collect();
    pts.extend(pts.first().copied());
    BSplineCurve::interpolate_points(&pts, None, 3, None, None, true, ParamMethod::ChordLength)
        .unwrap()
}

fn nurb(c: &BSplineCurve) -> Nurb {
    let g = c.gw_nurb().expect("section without NURB");
    Nurb {
        p: g.degree(),
        cpts: g.cpts().iter().map(|q| [q.x, q.y, q.z, q.w]).collect(),
        knots: g.knots().to_vec(),
    }
}

fn skin(from: usize, to: usize) -> Skin {
    Skin {
        sections: (from..to).map(|i| nurb(&section(i))).collect(),
        degree: 1,
        cap: true,
    }
}

fn rust_skin(k: &mut Kernel, from: usize, to: usize) -> Brep {
    let cs: Vec<Curve> = (from..to).map(|i| Curve::from(section(i))).collect();
    let refs: Vec<Option<&Curve>> = cs.iter().map(Some).collect();
    create_skin_primitive(k, &vec![1; cs.len()], &refs, 1, true)
        .unwrap()
        .brep
}

fn rust_shape(r: solidcore::Result<Brep>) -> String {
    match r {
        Ok(b) => format!(
            "F/E/V {}/{}/{}",
            b.face_list().len(),
            b.edge_list().len(),
            b.vertex_list().len()
        ),
        Err(e) => format!("{e:?}"),
    }
}

#[test]
#[ignore = "one-off #1036 diagnostic: prints the C++ (v1.0.1) and Rust results of the three Booleans"]
fn issue_1036_cpp_booleans() {
    let (a, b) = (skin(0, 5), skin(2, 7));
    let mut report = Vec::new();
    let measured = ffi::run(Op::Measure, &a, &b).unwrap();
    report.push(format!(
        "C++ skins A, B [status F E V vstatus volume]: {measured:?}"
    ));
    for (name, op) in [
        ("A - B", Op::Difference),
        ("A n B", Op::Intersection),
        ("A u B", Op::Union),
    ] {
        let out = ffi::run(op, &a, &b).unwrap();
        // the two skin records come first: 6 values each on success
        report.push(format!(
            "C++ {name} [skins.., status has F E V vstatus volume]: {out:?}"
        ));
        let mut k = Kernel::default();
        let (ra, rb) = (rust_skin(&mut k, 0, 5), rust_skin(&mut k, 2, 7));
        let r = match op {
            Op::Difference => boolean_difference(&mut k, ra, rb),
            Op::Intersection => boolean_intersection(&mut k, ra, rb),
            _ => boolean_union(&mut k, ra, rb),
        };
        report.push(format!("Rust {name}: {}", rust_shape(r)));
    }
    for line in &report {
        println!("{line}");
    }
    assert!(!measured.is_empty());
}
