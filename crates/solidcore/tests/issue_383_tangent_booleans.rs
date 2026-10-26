//! tangency-labs/solidcore#383: a cylinder or torus tangent to a planar face along a line
//! or circle. Upstream's `SmApiBoolean` returns a body that is not a closed solid, and
//! exact mass properties then fail on it. The cases are the issue's table (box and a
//! cylinder or torus, millimetres).
use solidcore::boolean::api_brep::{boolean_difference, boolean_intersection, boolean_union};
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::math::{Point3, Vec3};
use solidcore::primitive::api_primitives::{create_box, create_cylinder, create_torus};
use solidcore::topology::brep::p1;
use solidcore::topology::brep::p2::{is_manifold_solid, shell_vertices};
use solidcore::topology::edge;
use solidcore::topology::proto_topology::p1::api_transform::{
    api_rotate, api_translate, TransformObject,
};
use solidcore::{Brep, Kernel};

fn volume(k: &Kernel, b: &mut Brep) -> solidcore::Result<f64> {
    let store = AttributeStore::default();
    let attrs = PropertyAttributes {
        store: &store,
        key_of: &|_| ObjectKey(0),
    };
    brep_volume(k, b, 1e-6, &attrs)
}

/// The box with corners `a` and `b`.
fn cuboid(k: &mut Kernel, a: [f64; 3], b: [f64; 3]) -> solidcore::Result<Brep> {
    create_box(
        k,
        Point3::new(a[0], a[1], a[2]),
        b[0] - a[0],
        b[1] - a[1],
        b[2] - a[2],
    )
}

/// A cylinder of radius `r` and height `h` from `base` along the axis `axis` (0 X,
/// 1 Y, 2 Z).
fn cylinder(k: &mut Kernel, axis: usize, base: [f64; 3], r: f64, h: f64) -> solidcore::Result<Brep> {
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

fn torus(k: &mut Kernel, c: [f64; 3], big: f64, small: f64) -> solidcore::Result<Brep> {
    create_torus(k, Point3::new(c[0], c[1], c[2]), big, small)
}

/// What `is_manifold_solid` found wrong: the region count, the shell vertices and the
/// edges that are not manifold (with their edge use counts).
fn why_not_manifold(b: &Brep) -> String {
    let regions = p1::regions(b, None).map(|r| r.len());
    let shell_vertices = shell_vertices(b, None).map(|v| v.len());
    let bad: Vec<usize> = p1::edges(b, None)
        .unwrap_or_default()
        .into_iter()
        .filter(|&e| !edge::is_manifold(b, e).unwrap_or(false))
        .map(|e| b.edge(e).map(|r| r.edgeuses.len()).unwrap_or(0))
        .collect();
    format!(
        "not a manifold solid: regions {regions:?}, shell vertices {shell_vertices:?}, \
         {} non-manifold edges (edge uses {bad:?}) of {:?}",
        bad.len(),
        p1::edges(b, None).map(|e| e.len())
    )
}

#[derive(Clone, Copy, Debug)]
enum Op {
    Union,
    Cut,
    Intersect,
}

fn boolean(k: &mut Kernel, op: Op, x: Brep, y: Brep) -> Result<Brep, String> {
    match op {
        Op::Union => boolean_union(k, x, y),
        Op::Cut => boolean_difference(k, x, y),
        Op::Intersect => boolean_intersection(k, x, y),
    }
    .map_err(|e| format!("{op:?}: {e:?}"))
}

/// `op` of `a` and `b` (with the #746 slit-loop removal when `slits`), checked by volume
/// against inclusion-exclusion (`|A| + |B| - |A n B|`, `|A| - |A n B|`, `|A n B|`, the
/// intersection's volume from its own Boolean); the result's manifold status is
/// reported, not required (a cut tangent from inside leaves a knife edge).
fn run(
    op: Op,
    a: impl Fn(&mut Kernel) -> solidcore::Result<Brep>,
    b: impl Fn(&mut Kernel) -> solidcore::Result<Brep>,
    slits: bool,
) -> Result<String, String> {
    let e = |x: solidcore::SolidError| format!("{x:?}");
    let mut k = Kernel::default();
    k.fixes.remove_slit_loops_where_solids_touch = slits;
    let (x, y) = (a(&mut k).map_err(e)?, b(&mut k).map_err(e)?);
    let mut out = boolean(&mut k, op, x, y)?;
    let got = volume(&k, &mut out).map_err(|x| format!("volume: {x:?}"))?;
    let manifold = is_manifold_solid(&out).map_err(e)?;
    let shape = if manifold {
        String::from("manifold")
    } else {
        why_not_manifold(&out)
    };
    let (x, y) = (a(&mut k).map_err(e)?, b(&mut k).map_err(e)?);
    let mut both = boolean(&mut k, Op::Intersect, x, y)?;
    let vi = volume(&k, &mut both).map_err(|x| format!("intersection volume: {x:?}"))?;
    let (mut x, mut y) = (a(&mut k).map_err(e)?, b(&mut k).map_err(e)?);
    let va = volume(&k, &mut x).map_err(e)?;
    let vb = volume(&k, &mut y).map_err(e)?;
    let want = match op {
        Op::Union => va + vb - vi,
        Op::Cut => va - vi,
        Op::Intersect => vi,
    };
    if got.is_finite() && (got - want).abs() <= 1e-5 * want.abs().max(1.0) {
        Ok(format!("V {got} ({shape})"))
    } else {
        Err(format!("V {got}, want {want} ({shape})"))
    }
}

#[test]
fn issue_383_tangent_cylinder_or_torus_gives_a_closed_solid() {
    let mut report = Vec::new();
    for slits in [false, true] {
    let cases: Vec<(&str, Result<String, String>)> = vec![
        ("seed 1", run(Op::Union, |k| cylinder(k, 2, [15.0, 25.0, 22.5], 7.5, 22.5), |k| cuboid(k, [22.5, 22.5, 30.0], [27.5, 50.0, 45.0]), slits)),
        ("seed 4", run(Op::Cut, |k| cuboid(k, [22.5, 25.0, 22.5], [30.0, 32.5, 50.0]), |k| torus(k, [20.0, 22.5, 25.0], 10.0, 2.5), slits)),
        ("seed 32", run(Op::Cut, |k| cuboid(k, [7.5, 15.0, 0.0], [25.0, 32.5, 27.5]), |k| cylinder(k, 2, [12.5, 20.0, 10.0], 5.0, 30.0), slits)),
        ("seed 36", run(Op::Union, |k| cylinder(k, 1, [5.0, 30.0, 30.0], 2.5, 7.5), |k| cuboid(k, [2.5, 17.5, 20.0], [27.5, 42.5, 27.5]), slits)),
        ("seed 39", run(Op::Cut, |k| cuboid(k, [20.0, 5.0, 5.0], [35.0, 20.0, 32.5]), |k| cylinder(k, 2, [25.0, 12.5, 22.5], 5.0, 25.0), slits)),
        ("seed 65", run(Op::Union, |k| cuboid(k, [12.5, 2.5, 30.0], [37.5, 32.5, 55.0]), |k| torus(k, [25.0, 7.5, 25.0], 10.0, 5.0), slits)),
        ("seed 74 intersect", run(Op::Intersect, |k| cuboid(k, [7.5, 12.5, 12.5], [17.5, 32.5, 22.5]), |k| cylinder(k, 1, [10.0, 0.0, 22.5], 10.0, 27.5), slits)),
        ("seed 74 cut", run(Op::Cut, |k| cuboid(k, [7.5, 12.5, 12.5], [17.5, 32.5, 22.5]), |k| cylinder(k, 1, [10.0, 0.0, 22.5], 10.0, 27.5), slits)),
        ("seed 85 intersect", run(Op::Intersect, |k| cuboid(k, [17.5, 5.0, 7.5], [32.5, 25.0, 37.5]), |k| cylinder(k, 0, [15.0, 10.0, 30.0], 5.0, 25.0), slits)),
        ("seed 85 cut", run(Op::Cut, |k| cuboid(k, [17.5, 5.0, 7.5], [32.5, 25.0, 37.5]), |k| cylinder(k, 0, [15.0, 10.0, 30.0], 5.0, 25.0), slits)),
        ("seed 86", run(Op::Cut, |k| cuboid(k, [30.0, 20.0, 15.0], [35.0, 30.0, 32.5]), |k| torus(k, [20.0, 17.5, 30.0], 12.5, 2.5), slits)),
    ];
    report.push((slits, cases));
    }
    let failed = report
        .iter()
        .flat_map(|(slits, cases)| cases.iter().map(move |c| (slits, c)))
        .filter(|(_, (_, r))| r.is_err())
        .count();
    assert_eq!(failed, 0, "{report:#?}");
}

