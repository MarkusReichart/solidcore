//! Probe (not a regression test) for `corpus.parity`'s volume miss on
//! `Fillet/sphere_diff_3cylinders.smb`: per face, Rust next to C++, the area as
//! `corpus.parity` integrates it and the signed area and volume of
//! `ComputePreciseProperties(SM_OT_SAME, ..)` with the whole B-rep's origin and area
//! estimate (as `SmBrep::ComputePreciseProperties`), so the face or faces carrying the
//! whole-B-rep difference show. Run with
//! `cargo test -p difftest --features cpp --test probe_sphere_diff_faces -- --ignored --nocapture`
//! (`PROBE_FILE` names another model under `upstream/TestFiles/pt_TestFiles`,
//! `PROBE_ACCURACY` another accuracy). It asserts nothing.
#![cfg(feature = "cpp")]
#![allow(
    clippy::unwrap_used,
    clippy::indexing_slicing,
    clippy::print_stdout,
    clippy::panic
)]

use std::collections::BTreeMap;

use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::object::ObjectRef;
use solidcore::kernel::Kernel;
use solidcore::mass::integrator::properties::face_properties;
use solidcore::nurbs::b_spline_surface_n_lib::p2::OrientType;
use solidcore::query::brep::BoundingBoxOptions;
use solidcore::surface::surface::SurfaceKind;
use solidcore::topology::brep::p1 as bl;
use solidcore::topology::brep::p10::TopologyEnv;
use solidcore::topology::brep::p7::read_from_file;
use solidcore::topology::database_io::FileType;
use solidcore::topology::face::p2::{PPF_ALL, PPF_AREA};
use solidcore::curve::curve::CurveKind;
use solidcore::topology::face::p1::create_oriented_trim_boundary;
use solidcore::topology::Brep;
use solidcore::EntityRef;
use solidcore_ffi::units::corpus_faces as ffi;

const CORPUS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../upstream/TestFiles/pt_TestFiles");

fn key_of(e: EntityRef) -> ObjectRef {
    let (type_code, kind, i) = match e {
        EntityRef::Brep(x) => (16011, 0u64, x.index()),
        EntityRef::Region(x) => (16012, 1, x.index()),
        EntityRef::Shell(x) => (16013, 2, x.index()),
        EntityRef::Faceuse(x) => (16022, 3, x.index()),
        EntityRef::Face(x) => (16002, 4, x.index()),
        EntityRef::Loopuse(x) => (16004, 5, x.index()),
        EntityRef::Loop(x) => (16003, 6, x.index()),
        EntityRef::Edgeuse(x) => (16005, 7, x.index()),
        EntityRef::Edge(x) => (16016, 8, x.index()),
        EntityRef::Vertexuse(x) => (16006, 9, x.index()),
        EntityRef::Vertex(x) => (16021, 10, x.index()),
    };
    ObjectRef {
        type_code,
        key: ObjectKey((1 << 36) | ((kind + 1) << 32) | u64::from(i)),
    }
}

fn kind_name(k: &SurfaceKind) -> &'static str {
    match k {
        SurfaceKind::Plane(_) => "plane",
        SurfaceKind::Cone(_) => "cone",
        SurfaceKind::Sphere(_) => "sphere",
        SurfaceKind::Torus(_) => "torus",
        SurfaceKind::BSpline(_) => "bspline",
        SurfaceKind::SurfOfRevolution(_) => "revolution",
        SurfaceKind::SurfOfExtrusion(_) => "extrusion",
        SurfaceKind::OffsetSurface(_) => "offset",
        _ => "other",
    }
}

fn read(k: &mut Kernel, path: &str) -> Brep {
    let bytes = std::fs::read(path).unwrap();
    let mut b = Brep::new_with_infinite_region().unwrap();
    let mut attributes = AttributeStore::default();
    let mut records = BTreeMap::new();
    let mut env = TopologyEnv {
        all_attributes: &[],
        attributes: &mut attributes,
        key_of: &key_of,
        records: &mut records,
        next_key: 1 << 62,
    };
    read_from_file(k, &mut b, bytes, FileType::Ascii, false, &mut env).unwrap();
    b
}

/// Per face: (kind, area, signed area, signed volume, face box midpoint).
fn rust_faces(path: &str, accuracy: f64) -> Vec<(String, f64, f64, f64, [f64; 3])> {
    let mut k = Kernel::default();
    let mut b = read(&mut k, path);
    let bbox =
        solidcore::query::brep::calculate_bounding_box(&k, &mut b, BoundingBoxOptions::default())
            .unwrap();
    let s = bbox.size();
    let est = s.z.max(s.x.max(s.y));
    let est = est * est / 100.0;
    let mut out = Vec::new();
    for f in bl::faces(&b, None).unwrap() {
        let fbox = solidcore::query::face::p1::calculate_bounding_box(&k, &b, f, false).unwrap();
        let fs = fbox.size();
        let d = fs.x.max(fs.y).max(fs.z);
        let area = face_properties(
            &k,
            &mut b,
            f,
            OrientType::Same,
            accuracy,
            fbox.mid(),
            d * d / 100.0,
            0.0,
            PPF_AREA,
        )
        .unwrap()
        .area
        .abs();
        let all = face_properties(
            &k,
            &mut b,
            f,
            OrientType::Same,
            accuracy,
            bbox.mid(),
            est,
            0.0,
            PPF_ALL,
        )
        .unwrap();
        let kind = b
            .face(f)
            .unwrap()
            .surface
            .as_ref()
            .map_or("none", |s| kind_name(s.kind()));
        let m = fbox.mid();
        out.push((kind.to_string(), area, all.area, all.volume, [m.x, m.y, m.z]));
    }
    out
}

#[test]
#[ignore = "probe for corpus.parity's sphere_diff_3cylinders volume miss; run on request"]
fn probe_sphere_diff_faces() {
    let file = std::env::var("PROBE_FILE")
        .unwrap_or_else(|_| "Fillet/sphere_diff_3cylinders.smb".to_string());
    let accuracy: f64 = std::env::var("PROBE_ACCURACY")
        .ok()
        .and_then(|a| a.parse().ok())
        .unwrap_or(1e-9);
    let path = format!("{CORPUS}/{file}");
    let rust = rust_faces(&path, accuracy);
    let cpp = ffi::measure(&path, accuracy).unwrap();
    assert_eq!(cpp[0], 0.0, "C++ read status");
    let n = cpp[1] as usize;
    println!("{file} at accuracy {accuracy:e}: Rust {} faces, C++ {n}", rust.len());
    println!(
        "{:>3} {:>10} {:>8} {:>20} {:>20} {:>12} {:>20} {:>20} {:>12} {:>12} {:>9}",
        "i", "rust kind", "cpp type", "rust area", "cpp area", "d area", "rust vol", "cpp vol",
        "d vol", "d area_same", "mid gap"
    );
    let (mut rv, mut cv, mut ra, mut ca) = (0.0, 0.0, 0.0, 0.0);
    for (i, r) in rust.iter().enumerate() {
        let Some(c) = cpp.get(2 + i * ffi::PER_FACE..2 + (i + 1) * ffi::PER_FACE) else {
            println!("{i:>3} no C++ face");
            continue;
        };
        let gap = ((r.4[0] - c[4]).powi(2) + (r.4[1] - c[5]).powi(2) + (r.4[2] - c[6]).powi(2))
            .sqrt();
        println!(
            "{i:>3} {:>10} {:>8} {:>20.12} {:>20.12} {:>12.3e} {:>20.12} {:>20.12} {:>12.3e} {:>12.3e} {:>9.2e}  mid ({:.3}, {:.3}, {:.3}) cpp status {}",
            r.0, c[0], r.1, c[1], r.1 - c[1], r.3, c[3], r.3 - c[3], r.2 - c[2], gap,
            r.4[0], r.4[1], r.4[2], c[7]
        );
        rv += r.3;
        cv += c[3];
        ra += r.1;
        ca += c[1];
    }
    println!("sum area: Rust {ra:.12} C++ {ca:.12} diff {:.6e}", ra - ca);
    println!("sum signed SAME volume: Rust {rv:.12} C++ {cv:.12} diff {:.6e}", rv - cv);
}

/// The trim boundary of face `face` the integrator uses, as the shim's layout:
/// `loops`, then per curve `(degree, t0, t1, control points (x, y, w))`.
type Trim = (usize, f64, f64, Vec<[f64; 3]>);

fn rust_trims(path: &str, face: usize) -> (Vec<usize>, Vec<Trim>) {
    let mut k = Kernel::default();
    let mut b = read(&mut k, path);
    let f = bl::faces(&b, None).unwrap()[face];
    let key = |e: EntityRef| ObjectRef {
        type_code: 0,
        key: ObjectKey(match e {
            EntityRef::Edgeuse(eu) => u64::from(eu.index()),
            _ => 0,
        }),
    };
    let t = create_oriented_trim_boundary(&k, &mut b, f, true, true, None, &key).unwrap();
    let curves = t
        .curves
        .iter()
        .map(|c| {
            // upstream's trims are SmBSplineCurve*: the analytic kinds derive from it, so
            // their B-spline base is what C++ hands back
            // (a lazy line makes its NURB on a copy)
            let (name, base): (&str, Option<std::borrow::Cow<'_, _>>) = match c.kind() {
                CurveKind::BSpline(b) => ("bspline", Some(std::borrow::Cow::Borrowed(b))),
                CurveKind::Line(l) => ("line", Some(l.base_with_nurb().unwrap())),
                CurveKind::Ellipse(e) => ("ellipse", Some(std::borrow::Cow::Borrowed(e.base()))),
                CurveKind::Hyperbola(h) => {
                    ("hyperbola", Some(std::borrow::Cow::Borrowed(h.base())))
                }
                CurveKind::Parabola(p) => ("parabola", Some(std::borrow::Cow::Borrowed(p.base()))),
                _ => ("other", None),
            };
            let Some(bs) = base else {
                println!("  Rust trim of kind {name} has no B-spline base: {c:?}");
                return (0, f64::NAN, f64::NAN, Vec::new());
            };
            if name != "bspline" {
                println!("  Rust trim kind {name}");
            }
            let (cps, w) = bs.control_polygon().unwrap();
            let ivl = bs.natural_interval();
            let pts = cps
                .iter()
                .enumerate()
                .map(|(j, p)| [p.x, p.y, w.get(j).copied().unwrap_or(1.0)])
                .collect();
            (bs.degree().unwrap(), ivl.min(), ivl.max(), pts)
        })
        .collect();
    (t.loops, curves)
}

fn cpp_trims(path: &str, face: usize) -> (Vec<usize>, Vec<Trim>) {
    let v = ffi::trims(path, face as i64).unwrap();
    // SM_SUCCESS is 1000
    assert_eq!(v[0], 1000.0, "C++ trim boundary status");
    let nl = v[1] as usize;
    let loops: Vec<usize> = v[2..2 + nl].iter().map(|&x| x as usize).collect();
    let mut i = 2 + nl;
    let mut curves = Vec::new();
    while i < v.len() {
        let (deg, n, t0, t1) = (v[i] as usize, v[i + 1] as usize, v[i + 2], v[i + 3]);
        i += 4;
        let pts = (0..n).map(|j| [v[i + 3 * j], v[i + 3 * j + 1], v[i + 3 * j + 2]]).collect();
        i += 3 * n;
        curves.push((deg, t0, t1, pts));
    }
    (loops, curves)
}

#[test]
#[ignore = "probe for corpus.parity's sphere_diff_3cylinders volume miss; run on request"]
fn probe_sphere_diff_trims() {
    let file = std::env::var("PROBE_FILE")
        .unwrap_or_else(|_| "Fillet/sphere_diff_3cylinders.smb".to_string());
    let faces: Vec<usize> = std::env::var("PROBE_FACES")
        .unwrap_or_else(|_| "7,6".to_string())
        .split(',')
        .map(|s| s.trim().parse().unwrap())
        .collect();
    let path = format!("{CORPUS}/{file}");
    for face in faces {
        let (rl, rc) = rust_trims(&path, face);
        let (cl, cc) = cpp_trims(&path, face);
        println!("face {face}: loops Rust {rl:?} C++ {cl:?}; curves Rust {} C++ {}", rc.len(), cc.len());
        for (j, (r, c)) in rc.iter().zip(&cc).enumerate() {
            let same_shape = r.0 == c.0 && r.3.len() == c.3.len();
            let mut max = 0.0_f64;
            if same_shape {
                for (p, q) in r.3.iter().zip(&c.3) {
                    for d in 0..3 {
                        max = max.max((p[d] - q[d]).abs());
                    }
                }
            }
            println!(
                "  curve {j}: Rust deg {} n {} [{:.15}, {:.15}] | C++ deg {} n {} [{:.15}, {:.15}] | max cp diff {}",
                r.0, r.3.len(), r.1, r.2, c.0, c.3.len(), c.1, c.2,
                if same_shape { format!("{max:.3e}") } else { "shape differs".to_string() }
            );
            if !same_shape || max > 0.0 {
                println!("    Rust first/last cp {:?} {:?}", r.3.first(), r.3.last());
                println!("    C++  first/last cp {:?} {:?}", c.3.first(), c.3.last());
            }
        }
    }
}

/// One direct `CreateUVTrimCurve` result: `(status, degree, n, t0, t1, e0, e1, first,
/// last)`, status 0 for success on the Rust side and 1000 (`SM_SUCCESS`) on the C++.
type DropRes = (f64, usize, usize, f64, f64, f64, f64, [f64; 2], [f64; 2]);

fn rust_drops(path: &str, face: usize) -> Vec<DropRes> {
    use solidcore::topology::edgeuse::p1 as eu1;
    use solidcore::topology::edgeuse::p3::create_uv_trim_curve;
    let mut k = Kernel::default();
    let mut b = read(&mut k, path);
    let f = bl::faces(&b, None).unwrap()[face];
    let mut out = Vec::new();
    for l in b.face(f).unwrap().loops.clone() {
        let lu = b.loop_(l).unwrap().loopuse.unwrap();
        for eu in eu1::loopuse_edgeuses(&b, lu).unwrap() {
            let e = eu1::edge(&b, eu).unwrap().unwrap();
            let eivl = b.edge(e).unwrap().interval;
            match create_uv_trim_curve(&k, &mut b, eu, false, None) {
                Ok(made) => {
                    let base = match made.curve.kind() {
                        CurveKind::BSpline(c) => Some(std::borrow::Cow::Borrowed(c)),
                        CurveKind::Line(l) => Some(l.base_with_nurb().unwrap()),
                        _ => None,
                    };
                    let Some(bs) = base else {
                        println!("  Rust drop of edge use {eu:?}: not a B-spline: {:?}", made.curve);
                        continue;
                    };
                    let (cps, _) = bs.control_polygon().unwrap();
                    let ivl = bs.natural_interval();
                    let (p0, pn) = (cps[0], cps[cps.len() - 1]);
                    out.push((
                        0.0,
                        bs.degree().unwrap(),
                        cps.len(),
                        ivl.min(),
                        ivl.max(),
                        eivl.min(),
                        eivl.max(),
                        [p0.x, p0.y],
                        [pn.x, pn.y],
                    ));
                }
                Err(err) => {
                    println!("  Rust drop of edge use {eu:?} failed: {err}");
                    out.push((1.0, 0, 0, 0.0, 0.0, eivl.min(), eivl.max(), [0.0; 2], [0.0; 2]));
                }
            }
        }
    }
    out
}

fn cpp_drops(path: &str, face: usize) -> Vec<DropRes> {
    let v = ffi::drops(path, face as i64).unwrap();
    assert_eq!(v[0], 0.0, "C++ read status");
    v[1..]
        .chunks(11)
        .map(|c| {
            (
                c[0],
                c[1] as usize,
                c[2] as usize,
                c[3],
                c[4],
                c[5],
                c[6],
                [c[7], c[8]],
                [c[9], c[10]],
            )
        })
        .collect()
}

/// `SmEdgeuse::CreateUVTrimCurve` called directly on every edge use of the face, on both
/// sides, before `sm_FixUVTrimCurves` touches anything: whether the single-span trim
/// of face 7's edge 30 comes from the drop itself or from the corner fix after it.
#[test]
#[ignore = "probe for corpus.parity's sphere_diff_3cylinders volume miss; run on request"]
fn probe_sphere_diff_drops() {
    let file = std::env::var("PROBE_FILE")
        .unwrap_or_else(|_| "Fillet/sphere_diff_3cylinders.smb".to_string());
    let faces: Vec<usize> = std::env::var("PROBE_FACES")
        .unwrap_or_else(|_| "7,6".to_string())
        .split(',')
        .map(|s| s.trim().parse().unwrap())
        .collect();
    let path = format!("{CORPUS}/{file}");
    for face in faces {
        let r = rust_drops(&path, face);
        let c = cpp_drops(&path, face);
        println!("face {face}: edge uses Rust {} C++ {}", r.len(), c.len());
        for (i, (r, c)) in r.iter().zip(&c).enumerate() {
            println!(
                "  eu {i}: edge [{:.12}, {:.12}] | Rust st {} deg {} n {} [{:.12}, {:.12}] {:?}->{:?} | C++ st {} deg {} n {} [{:.12}, {:.12}] {:?}->{:?}",
                r.5, r.6, r.0, r.1, r.2, r.3, r.4, r.7, r.8, c.0, c.1, c.2, c.3, c.4, c.7, c.8
            );
        }
    }
}

/// The control point count of edge use `eu`'s UV trim curve (`-1` without one).
fn trim_points(b: &Brep, eu: solidcore::topology::EdgeuseId) -> f64 {
    use solidcore::topology::edgeuse::p1 as eu1;
    match eu1::uv_trim_curve(b, eu).unwrap().0 {
        None => -1.0,
        Some(c) => match c.kind() {
            CurveKind::BSpline(bs) => bs.control_polygon().unwrap().0.len() as f64,
            CurveKind::Line(l) => l.base_with_nurb().unwrap().control_polygon().unwrap().0.len() as f64,
            _ => -2.0,
        },
    }
}

/// The Rust side of `sc_corpus_faces_fixsteps`, in its layout.
fn rust_fixsteps(path: &str, face: usize) -> Vec<f64> {
    use solidcore::curve::curve::Curve;
    use solidcore::geom::global_solver::SolverOperation;
    use solidcore::math::{Point3, EFF_ZERO};
    use solidcore::topology::edgeuse::p1 as eu1;
    use solidcore::topology::edgeuse::p3::{get_or_create_uv_trim_curve, normalized_evaluate};
    use solidcore::topology::face::p2::fix_uv_trim_curves;
    use solidcore::topology::vertexuse as vu1;
    let mut k = Kernel::default();
    let mut b = read(&mut k, path);
    let f = bl::faces(&b, None).unwrap()[face];
    let uv_domain = b.face(f).unwrap().uv_domain;
    // GetUpwardLoopuses, as CreateUVTrimCurves takes them
    let lus = solidcore::topology::face::p4::upward_loopuses(&b, f, false, None).unwrap();
    let mut out = vec![0.0, f64::NAN];
    for &lu in &lus {
        for eu in eu1::loopuse_edgeuses(&b, lu).unwrap() {
            if eu1::edge(&b, eu).unwrap().is_none() || eu1::uv_trim_curve(&b, eu).unwrap().0.is_some() {
                continue;
            }
            let _ = get_or_create_uv_trim_curve(&k, &mut b, eu, true, false, true).unwrap();
        }
    }
    for (l, &lu) in lus.iter().enumerate() {
        for (j, eu) in eu1::loopuse_edgeuses(&b, lu).unwrap().into_iter().enumerate() {
            let cw = eu1::cw_edgeuse(&b, eu).unwrap().unwrap();
            let vu = eu1::vertexuse(&b, eu).unwrap().unwrap();
            let v = vu1::vertex(&b, vu).unwrap().unwrap();
            let p0 = normalized_evaluate(&k, &mut b, eu, 0.0, true, false, true).unwrap().point;
            let p1 = normalized_evaluate(&k, &mut b, cw, 1.0, true, false, true).unwrap().point;
            let gap = p0.distance(p1);
            let scaled = EFF_ZERO * (1.0 + p0.max_dimension());
            let nan = f64::NAN;
            if !(gap > scaled) {
                out.extend([l as f64, j as f64, gap, 0.0, nan, nan, nan, nan, nan, nan,
                    trim_points(&b, eu), trim_points(&b, cw), nan, nan, nan, nan, nan, nan, nan,
                    nan]);
                continue;
            }
            let cur: Curve = eu1::uv_trim_curve(&b, eu).unwrap().0.unwrap().clone();
            let ivl = cur.natural_interval();
            let same = eu1::orientation(&b, eu).unwrap() == OrientType::Same;
            let t_curr = if same { ivl.min() } else { ivl.max() };
            let uv_curr = cur.evaluate_point(&k, t_curr).unwrap();
            let uv_vertex = match vu1::vertex_face_gap(&k, &mut b, vu, false).unwrap() {
                Some(g) => g.face_uv,
                None => vu1::compute_uv_point(&k, &mut b, vu, true).unwrap(),
            };
            let uv3 = Point3::new(uv_vertex.x, uv_vertex.y, 0.0);
            let edge_gap = uv3.distance(uv_curr);
            let d = cur
                .drop_point(&k, &ivl, uv3, None, 1.1 * (EFF_ZERO + edge_gap), Some(t_curr), SolverOperation::Intersect)
                .ok();
            let prev: Curve = eu1::uv_trim_curve(&b, cw).unwrap().0.unwrap().clone();
            let pivl = prev.natural_interval();
            let t_prev = if eu1::orientation(&b, cw).unwrap() == OrientType::Same {
                pivl.max()
            } else {
                pivl.min()
            };
            let uv_prev = prev.evaluate_point(&k, t_prev).unwrap();
            let prev_gap = uv3.distance(uv_prev);
            let pd = prev
                .drop_point(&k, &pivl, uv3, None, 1.1 * (EFF_ZERO + prev_gap), Some(t_prev), SolverOperation::Intersect)
                .ok();
            let (n_cur, n_cw) = (trim_points(&b, eu), trim_points(&b, cw));
            let fixed = fix_uv_trim_curves(&k, &mut b, eu, cw, v, 0.0, gap, &uv_domain).unwrap_or(false);
            out.extend([
                l as f64,
                j as f64,
                gap,
                1.0,
                uv_vertex.x,
                uv_vertex.y,
                t_curr,
                d.map_or(0.0, |d| f64::from(u8::from(d.success))),
                d.map_or(nan, |d| d.parameter),
                d.map_or(nan, |d| d.distance),
                n_cur,
                n_cw,
                trim_points(&b, eu),
                trim_points(&b, cw),
                f64::from(u8::from(fixed)),
                0.0,
                t_prev,
                pd.map_or(0.0, |d| f64::from(u8::from(d.success))),
                pd.map_or(nan, |d| d.parameter),
                pd.map_or(nan, |d| d.distance),
            ]);
        }
    }
    out
}

/// `SmFace::CreateUVTrimCurves` replayed corner by corner on both sides: where
/// `sm_FixUVTrimCurves` first makes a trim differ (WR-449: face 7's edge 30 ends as a
/// single cubic span in C++, 46 points in the port, from the same 49-point drop).
#[test]
#[ignore = "probe for corpus.parity's sphere_diff_3cylinders volume miss; run on request"]
fn probe_sphere_diff_fixsteps() {
    let file = std::env::var("PROBE_FILE")
        .unwrap_or_else(|_| "Fillet/sphere_diff_3cylinders.smb".to_string());
    let faces: Vec<usize> = std::env::var("PROBE_FACES")
        .unwrap_or_else(|_| "7,6".to_string())
        .split(',')
        .map(|s| s.trim().parse().unwrap())
        .collect();
    let path = format!("{CORPUS}/{file}");
    for face in faces {
        let r = rust_fixsteps(&path, face);
        let c = ffi::fixsteps(&path, face as i64).unwrap();
        println!("face {face}: C++ read {} poles {}", c[0], c[1]);
        let names = [
            "loop", "j", "gap", "called", "uv.x", "uv.y", "t_curr", "drop ok", "drop t",
            "drop dist", "n cur", "n cw", "n cur after", "n cw after", "fixed", "-", "t_prev",
            "prev ok", "prev t", "prev dist",
        ];
        for (rr, cc) in r[2..].chunks(20).zip(c[2..].chunks(20)) {
            let row = |v: &[f64]| {
                names
                    .iter()
                    .zip(v)
                    .map(|(n, x)| format!("{n} {x:.15}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            };
            println!("  Rust: {}", row(rr));
            println!("  C++ : {}", row(cc));
        }
    }
}

/// Unique knots with multiplicities of a full knot vector (within 1e-15 relative).
fn unique_knots(full: &[f64]) -> Vec<(f64, usize)> {
    let mut out: Vec<(f64, usize)> = Vec::new();
    for &t in full {
        match out.last_mut() {
            Some((u, m)) if (t - *u).abs() <= 1e-15 * (1.0 + u.abs()) => *m += 1,
            _ => out.push((t, 1)),
        }
    }
    out
}

/// `(n, unique knots)` of a B-spline-based curve (a lazy line makes its NURB on a copy).
fn rust_shape(c: &solidcore::curve::curve::Curve) -> (usize, Vec<(f64, usize)>) {
    let base = match c.kind() {
        CurveKind::BSpline(b) => std::borrow::Cow::Borrowed(b),
        CurveKind::Line(l) => l.base_with_nurb().unwrap(),
        _ => panic!("not a B-spline"),
    };
    let n = base.control_polygon().unwrap().0.len();
    let (u, m) = base.knots(None).unwrap();
    (n, u.into_iter().zip(m).collect())
}

/// `sm_FixUVTrimCurves`' trim step alone, both sides, on the freshly made trim of edge use
/// `PROBE_EU` (default 2: face 7's edge 30) with the interval the fix trims it to at
/// face 7's corner 3 (`[min, tPrevNew]`, from `probe_sphere_diff_fixsteps`): the made
/// curve, `Trim` (with the interval it snaps to) and `EditParameterization` back.
#[test]
#[ignore = "probe for corpus.parity's sphere_diff_3cylinders volume miss; run on request"]
fn probe_sphere_diff_trimstep() {
    use solidcore::geom::extent1d::Extent1d;
    use solidcore::nurbs::b_spline_curve_n_lib::TrimMode;
    use solidcore::topology::edgeuse::p1 as eu1;
    use solidcore::topology::edgeuse::p3::create_uv_trim_curve;
    let env = |name: &str, default: f64| {
        std::env::var(name).ok().and_then(|v| v.parse().ok()).unwrap_or(default)
    };
    let file = std::env::var("PROBE_FILE")
        .unwrap_or_else(|_| "Fillet/sphere_diff_3cylinders.smb".to_string());
    let face = env("PROBE_FACE", 7.0) as usize;
    let eu_index = env("PROBE_EU", 2.0) as usize;
    let t0 = env("PROBE_T0", 0.499999957853150);
    let t1 = env("PROBE_T1", 0.749999985455706);
    let path = format!("{CORPUS}/{file}");

    // Rust
    let mut k = Kernel::default();
    let mut b = read(&mut k, &path);
    let f = bl::faces(&b, None).unwrap()[face];
    let mut all = Vec::new();
    for lu in solidcore::topology::face::p4::upward_loopuses(&b, f, false, None).unwrap() {
        all.extend(eu1::loopuse_edgeuses(&b, lu).unwrap());
    }
    let mut c = create_uv_trim_curve(&k, &mut b, all[eu_index], false, None).unwrap().curve;
    let made = rust_shape(&c);
    let full = c.natural_interval();
    let mut ivl = Extent1d::new();
    let _ = ivl.set_min_max(t0, t1);
    let trim = c.trim(&k, &mut ivl, TrimMode::Notify);
    let after_trim = rust_shape(&c);
    let edit = c.edit_parameterization(&k, &full);
    let after_edit = rust_shape(&c);

    // C++
    let v = ffi::trimstep(&path, face as i64, eu_index as i64, t0, t1).unwrap();
    assert_eq!(v[0], 0.0, "C++ status");
    let mut i = 1;
    let shape = |v: &[f64], i: &mut usize| {
        let (n, nk) = (v[*i] as usize, v[*i + 1] as usize);
        let knots = unique_knots(&v[*i + 2..*i + 2 + nk]);
        *i += 2 + nk;
        (n, knots)
    };
    let c_made = shape(&v, &mut i);
    let (c_trim_st, c_t0, c_t1) = (v[i], v[i + 1], v[i + 2]);
    i += 3;
    let c_after_trim = shape(&v, &mut i);
    let c_edit_st = v[i];
    i += 1;
    let c_after_edit = shape(&v, &mut i);

    println!("{file} face {face} edge use {eu_index}: trim to [{t0:.15}, {t1:.15}]");
    println!("made        Rust n {} knots {:?}", made.0, made.1);
    println!("made        C++  n {} knots {:?}", c_made.0, c_made.1);
    println!(
        "Trim        Rust {:?} -> [{:.15}, {:.15}] n {} knots {:?}",
        trim.as_ref().err(),
        ivl.min(),
        ivl.max(),
        after_trim.0,
        after_trim.1
    );
    println!(
        "Trim        C++  st {c_trim_st} -> [{c_t0:.15}, {c_t1:.15}] n {} knots {:?}",
        c_after_trim.0, c_after_trim.1
    );
    println!("EditParam   Rust {:?} n {} knots {:?}", edit.as_ref().err(), after_edit.0, after_edit.1);
    println!("EditParam   C++  st {c_edit_st} n {} knots {:?}", c_after_edit.0, c_after_edit.1);
}
