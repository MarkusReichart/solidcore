//! Probe (not a regression test) for `corpus.parity`'s volume miss on
//! `Fillet/sphere_diff_3cylinders.smb`: an estimate of each face's area and volume
//! contribution independent of the mass integrator, from the triangles of fine meshes
//! (the divergence theorem: a triangle (a, b, c) adds `(a-o)·((b-o)×(c-o))/6`, with `o`
//! the B-rep box midpoint `SmBrep::ComputePreciseProperties` uses), at three chord
//! tolerances, extrapolated to zero chord (an inscribed mesh's area and volume
//! deficits are linear in the chord tolerance, so `2 x(h/2) - x(h)` cancels the leading
//! term), next to the integrator's value for the same face. Run with
//! `cargo test --release -p solidcore-mesh --test probe_sphere_diff_mesh -- --ignored --nocapture`
//! (`PROBE_FILE` names another model under `upstream/TestFiles/pt_TestFiles`). It
//! asserts nothing.

#![allow(clippy::unwrap_used, clippy::indexing_slicing, clippy::print_stdout)]

use std::collections::BTreeMap;

use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::geom::object::ObjectRef;
use solidcore::kernel::Kernel;
use solidcore::mass::integrator::properties::face_properties;
use solidcore::math::Point3;
use solidcore::query::brep::BoundingBoxOptions;
use solidcore::topology::brep::p1 as bl;
use solidcore::topology::brep::p10::TopologyEnv;
use solidcore::topology::brep::p7::read_from_file;
use solidcore::topology::database_io::FileType;
use solidcore::topology::face::p2::PPF_ALL;
use solidcore::topology::{Brep, FaceId};
use solidcore::EntityRef;
use solidcore_mesh::{tessellate, OutputOptions, Quality, TessellateOptions, Tolerance};

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
        key: ObjectKey(((kind + 1) << 32) | u64::from(i)),
    }
}

/// Per face: (area, volume) of the mesh at relative chord tolerance `chord`.
fn mesh_faces(k: &Kernel, b: &Brep, origin: Point3, chord: f64) -> BTreeMap<FaceId, (f64, f64)> {
    let mut opts = TessellateOptions::preset(Quality::Fine);
    opts.chord_tolerance = Tolerance::Relative(chord);
    opts.max_triangles_per_face = 20_000_000;
    opts.max_refine_passes = 64;
    opts.output = OutputOptions {
        normals: false,
        uvs: false,
        face_ids: true,
        edge_polylines: false,
    };
    let mesh = tessellate(k, b, &opts).unwrap();
    if !mesh.failures.is_empty() || !mesh.unconverged_faces.is_empty() {
        println!(
            "chord {chord:e}: {} failures, unconverged faces {:?}",
            mesh.failures.len(),
            mesh.unconverged_faces
        );
    }
    let mut out = BTreeMap::new();
    for (t, f) in mesh.triangles.iter().zip(&mesh.triangle_faces) {
        let [a, bb, c] = t.map(|i| mesh.positions[i as usize] - origin);
        let n = (bb - a).cross(c - a);
        let e = out.entry(*f).or_insert((0.0, 0.0));
        e.0 += 0.5 * n.length();
        e.1 += a.dot(bb.cross(c)) / 6.0;
    }
    println!("chord {chord:e}: {} triangles", mesh.triangles.len());
    out
}

#[test]
#[ignore = "probe for corpus.parity's sphere_diff_3cylinders volume miss; run on request"]
fn probe_sphere_diff_mesh() {
    let file = std::env::var("PROBE_FILE")
        .unwrap_or_else(|_| "Fillet/sphere_diff_3cylinders.smb".to_string());
    let mut k = Kernel::default();
    let bytes = std::fs::read(format!("{CORPUS}/{file}")).unwrap();
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
    read_from_file(&mut k, &mut b, bytes, FileType::Ascii, false, &mut env).unwrap();
    let bbox =
        solidcore::query::brep::calculate_bounding_box(&k, &mut b, BoundingBoxOptions::default())
            .unwrap();
    let origin = bbox.mid();
    let s = bbox.size();
    let est = s.z.max(s.x.max(s.y));
    let est = est * est / 100.0;
    let faces = bl::faces(&b, None).unwrap();

    let levels = [4e-6, 2e-6, 1e-6];
    let meshes: Vec<_> = levels.iter().map(|&h| mesh_faces(&k, &b, origin, h)).collect();

    println!(
        "{file}: origin ({:.6}, {:.6}, {:.6}); per face, mesh at {levels:?}, extrapolated \
         (2 x(h/2) - x(h) of the last two), integrator (accuracy 1e-9, faceuse SAME)",
        origin.x, origin.y, origin.z
    );
    let (mut ta, mut tv, mut ia, mut iv) = (0.0, 0.0, 0.0, 0.0);
    for (i, f) in faces.iter().enumerate() {
        let at = |m: &BTreeMap<FaceId, (f64, f64)>| m.get(f).copied().unwrap_or((0.0, 0.0));
        let [m0, m1, m2] = [at(&meshes[0]), at(&meshes[1]), at(&meshes[2])];
        let ext = (2.0 * m2.0 - m1.0, 2.0 * m2.1 - m1.1);
        let p = face_properties(
            &k,
            &mut b,
            *f,
            solidcore::nurbs::b_spline_surface_n_lib::p2::OrientType::Same,
            1e-9,
            origin,
            est,
            0.0,
            PPF_ALL,
        )
        .unwrap();
        println!(
            "{i:>3} area {:>16.6} {:>16.6} {:>16.6} ext {:>16.6} int {:>16.6} | vol {:>16.6} {:>16.6} {:>16.6} ext {:>16.6} int {:>16.6} (|ext|-|int| {:>10.4})",
            m0.0, m1.0, m2.0, ext.0, p.area.abs(), m0.1, m1.1, m2.1, ext.1, p.volume,
            ext.1.abs() - p.volume.abs()
        );
        ta += ext.0;
        tv += ext.1;
        ia += p.area.abs();
        iv += p.volume.abs() * ext.1.signum();
    }
    println!("total area: mesh extrapolated {ta:.6}, integrator {ia:.6}");
    println!("total volume: mesh extrapolated {tv:.6}, integrator (signed as the mesh) {iv:.6}");
}

/// The faces the mesher leaves out of this model, with the error and the loop report
/// (`face_loop_report`) of each, at the Normal preset and at a fine chord.
#[test]
#[ignore = "probe for the faces solidcore-mesh fails on sphere_diff_3cylinders; run on request"]
fn probe_sphere_diff_mesh_failures() {
    let file = std::env::var("PROBE_FILE")
        .unwrap_or_else(|_| "Fillet/sphere_diff_3cylinders.smb".to_string());
    let mut k = Kernel::default();
    let bytes = std::fs::read(format!("{CORPUS}/{file}")).unwrap();
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
    read_from_file(&mut k, &mut b, bytes, FileType::Ascii, false, &mut env).unwrap();
    let faces = bl::faces(&b, None).unwrap();
    for (name, opts) in [
        ("normal", TessellateOptions::default()),
        ("fine chord 1e-6", {
            let mut o = TessellateOptions::preset(Quality::Fine);
            o.chord_tolerance = Tolerance::Relative(1e-6);
            o.max_triangles_per_face = 20_000_000;
            o
        }),
    ] {
        let mesh = tessellate(&k, &b, &opts).unwrap();
        println!("{file} {name}: {} triangles, {} failures", mesh.triangles.len(), mesh.failures.len());
        let edges = solidcore_mesh::tessellate_edge_set(&k, &b, &opts).unwrap();
        // a face that meshes, for comparison (index 2 by default)
        let good: usize = std::env::var("PROBE_GOOD_FACE")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(2);
        if let Some(&g) = faces.get(good) {
            println!("  reference face index {good}:");
            match solidcore_mesh::face_loop_report(&k, &b, g, &edges, &opts) {
                Ok(r) => println!("{r}"),
                Err(e) => println!("  loop report failed: {e}"),
            }
        }
        for fail in &mesh.failures {
            let index = match fail.entity {
                EntityRef::Face(f) => faces.iter().position(|g| *g == f),
                _ => None,
            };
            println!("  failure {:?} (face index {index:?}): {}", fail.entity, fail.error);
            if let EntityRef::Face(f) = fail.entity {
                match solidcore_mesh::face_loop_report(&k, &b, f, &edges, &opts) {
                    Ok(r) => println!("{r}"),
                    Err(e) => println!("  loop report failed: {e}"),
                }
            }
        }
    }
}
