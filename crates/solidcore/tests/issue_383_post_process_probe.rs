//! Probe (not for commit) for the torus_guard cases of `boolean.api_brep` where C++
//! measures a result and Rust has none: what Rust's Boolean returns.
use solidcore::boolean::api_brep::boolean_with_options;
use solidcore::boolean::merge::MergeOptions;
use solidcore::geom::a_object::BooleanOperation;
use solidcore::mass::integrator::properties::face_area;
use solidcore::math::Point3;
use solidcore::primitive::api_primitives::{create_box, create_torus};
use solidcore::topology::brep::p1;
use solidcore::Kernel;

/// The torus_guard cases where C++ measures a result and Rust has none: what Rust's
/// Boolean (and then its volume and face areas) return.
#[test]
fn torus_guard_spindle_and_tiny_probe() {
    use solidcore::geom::attribute::{AttributeStore, ObjectKey};
    use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
    let mut report = Vec::new();
    let cases = [
        (
            "spindle n",
            [1.0, -1.5, -1.5],
            [6.0, 3.0, 3.0],
            2.0,
            3.0,
            1.0,
            BooleanOperation::Intersection,
        ),
        (
            "spindle -",
            [1.0, -1.5, -1.5],
            [6.0, 3.0, 3.0],
            2.0,
            3.0,
            1.0,
            BooleanOperation::Difference,
        ),
        (
            "horn n",
            [1.0, -1.5, -1.5],
            [6.0, 3.0, 3.0],
            2.5,
            2.5,
            1.0,
            BooleanOperation::Intersection,
        ),
        (
            "tiny n",
            [0.5, -0.5, -1.0],
            [1.0, 1.0, 2.0],
            1.0,
            0.3,
            1e-4,
            BooleanOperation::Intersection,
        ),
        (
            "tiny -",
            [0.5, -0.5, -1.0],
            [1.0, 1.0, 2.0],
            1.0,
            0.3,
            1e-4,
            BooleanOperation::Difference,
        ),
    ];
    for (name, o, sz, major, minor, s, op) in cases {
        let mut k = Kernel::default();
        let a = create_box(
            &mut k,
            Point3::new(o[0] * s, o[1] * s, o[2] * s),
            sz[0] * s,
            sz[1] * s,
            sz[2] * s,
        )
        .unwrap();
        let t = create_torus(&mut k, Point3::new(0.0, 0.0, 0.0), major * s, minor * s);
        let Ok(t) = t else {
            report.push(format!("{name}: torus {t:?}"));
            continue;
        };
        let options = MergeOptions::new(false, false, true, false, false);
        match boolean_with_options(&mut k, a, t, op, options) {
            Ok(mut out) => {
                let store = AttributeStore::default();
                let attrs = PropertyAttributes {
                    store: &store,
                    key_of: &|_| ObjectKey(0),
                };
                let v = brep_volume(&k, &mut out.result, 1e-6, &attrs);
                let mut areas = Vec::new();
                for f in p1::faces(&out.result, None).unwrap() {
                    areas.push(format!("{:?}", face_area(&k, &mut out.result, f, 1e-4)));
                }
                report.push(format!("{name}: Ok, volume {v:?}, areas {areas:?}"));
            }
            Err(e) => report.push(format!("{name}: {e:?}")),
        }
    }
    // spindle n by stage: imprint only, then the whole Boolean without the post-process
    for (name, imprinting, post) in [
        ("imprint only", true, false),
        ("no post-process", false, false),
    ] {
        let mut k = Kernel::default();
        let a = create_box(&mut k, Point3::new(1.0, -1.5, -1.5), 6.0, 3.0, 3.0).unwrap();
        let t = create_torus(&mut k, Point3::new(0.0, 0.0, 0.0), 2.0, 3.0).unwrap();
        let options = MergeOptions::new(false, imprinting, post, false, false);
        let r = boolean_with_options(&mut k, a, t, BooleanOperation::Intersection, options);
        report.push(format!(
            "spindle n, {name}: {:?}",
            r.map(|o| p1::faces(&o.result, None).map(|f| f.len()))
        ));
    }
    panic!("{report:#?}");
}

/// The spindle torus (major 2, minor 3) against a box only, intersection and difference,
/// for a debugger run at the failing imprint (gspindle2.gdb): Rust's error for each.
#[test]
fn spindle_imprint_args_probe() {
    let mut report = Vec::new();
    for op in [BooleanOperation::Intersection, BooleanOperation::Difference] {
        let mut k = Kernel::default();
        let a = create_box(&mut k, Point3::new(1.0, -1.5, -1.5), 6.0, 3.0, 3.0).unwrap();
        let t = create_torus(&mut k, Point3::new(0.0, 0.0, 0.0), 2.0, 3.0).unwrap();
        let options = MergeOptions::new(false, false, true, false, false);
        let r = boolean_with_options(&mut k, a, t, op, options);
        report.push(format!(
            "{op:?}: {:?}",
            r.map(|o| p1::faces(&o.result, None).map(|f| f.len()))
        ));
    }
    panic!("{report:#?}");
}
