//! Filleting a boolean result (SolidCraft #324 `round_2_clear_of_hole`: a 40 x 20 x 10
//! block less a drill of radius 3 through it, then a radius 2 round on the top back
//! edge). The fillet's merge stops at `MergeCurveOnSurfaces` with
//! `SM_ERR_NULL_POINTER` (1006) when a face's surface is owned under a name the result's
//! topology does not resolve to. Run with `--nocapture` for each face's surface owner.

use solidcore::blend::fillet_executive::api::fillet_edges;
use solidcore::boolean::api_brep::boolean_difference;
use solidcore::error::EntityRef;
use solidcore::kernel::Kernel;
use solidcore::math::Point3;
use solidcore::primitive::api_primitives::{create_box, create_cylinder};
use solidcore::topology::brep::p3::{topology_entity, topology_ref};
use solidcore::topology::{edge, Brep, EdgeId};

/// Each face's surface owner: this B-rep's name for the face, another name of this
/// B-rep, or another B-rep's (decoded serial, entity) name.
fn owner_table(b: &Brep) -> Vec<String> {
    let serial = b.notify_serial();
    let mut rows = Vec::new();
    for (id, face) in b.faces().iter() {
        let owner = face.surface.as_ref().and_then(|s| s.owner());
        let row = match (owner, serial) {
            (None, _) => "no owner".to_string(),
            (Some(o), Some(s)) if o == topology_ref(s, EntityRef::Face(id)) => "own".to_string(),
            (Some(o), Some(s)) => match topology_entity(s, o.key) {
                Some(e) => format!("this B-rep's {e:?}"),
                None => {
                    let other = (o.key.0 >> 36) & 0x3ff_ffff;
                    format!(
                        "serial {other}: {:?} (key {:#x})",
                        topology_entity(other, o.key),
                        o.key.0
                    )
                }
            },
            (Some(o), None) => format!("unserialled B-rep; owner key {:#x}", o.key.0),
        };
        rows.push(format!("{id:?}: {row}"));
    }
    rows
}

fn edge_through(k: &Kernel, b: &Brep, p: Point3) -> EdgeId {
    use solidcore::geom::global_solver::SolverOperation;
    let found: Vec<EdgeId> = b
        .edges()
        .ids()
        .filter(|&e| {
            let curve = edge::curve(b, e).unwrap().unwrap();
            let ivl = edge::interval(b, e).unwrap();
            let d = curve
                .drop_point(k, &ivl, p, None, 1e-6, None, SolverOperation::Minimize)
                .unwrap();
            d.success && d.distance <= 1e-6
        })
        .collect();
    assert_eq!(found.len(), 1, "edges through {p:?}: {found:?}");
    found[0]
}

#[test]
fn fillet_a_difference_result() {
    let mut k = Kernel::default();
    let block = create_box(&mut k, Point3::new(0.0, 0.0, 0.0), 40.0, 20.0, 10.0).unwrap();
    let drill = create_cylinder(&mut k, Point3::new(20.0, 14.5, -1.0), 3.0, 12.0).unwrap();
    let mut body = boolean_difference(&mut k, block, drill).unwrap();
    for row in owner_table(&body) {
        println!("{row}");
    }
    let e = edge_through(&k, &body, Point3::new(10.0, 20.0, 10.0));
    fillet_edges(&mut k, &mut body, &[e], 2.0, 1, 1, 1.0).unwrap();
    assert!(solidcore::topology::brep::p2::is_manifold_solid(&body).unwrap());
    // every link and the cached topology fields the fillet API refreshes on return
    solidcore::topology::integrity::check(&body).unwrap();
}
