//! The end of `SmFilletExecutive::DoFilleting` on a [`FilletSession`]: the global merge,
//! the stitch and clean-up after the fillet faces are merged into the target, and the
//! session's `DoFillet` that `SmApiFilletEdges` drives.
//!
//! Upstream: `source/SMLib/src/SmFilletExecutive.cpp` (`DoFillet`, `DoFilleting` from the
//! global merge on), `source/SM_API/src/SmApiBrep.cpp` (`SmApiFilletEdges`).
//!
//! | Rust                                   | Upstream                                              |
//! |----------------------------------------|-------------------------------------------------------|
//! | [`FilletSession::global_merge`]        | `DoFilleting`: `m_vTI.IntersectInsertRelate()` when `m_bDoGlobalMerge` |
//! | [`FilletSession::finish_fillet_merge`] | `DoFilleting` after `MergeBrep`: delete the fillet B-rep, `SmStitch::DoStitching` (5 tries), `RemoveTopoEdges`, `BooleanPostProcess(TRUE)`, the faces' `Notify(SM_NO_PRE_EDIT)`, `OrientTrimmedSurfaces` |
//! | [`FilletSession::do_fillet`]           | `SmFilletExecutive::DoFillet` (`CreateFilletCorners`, `DoFilleting`) |
//! | [`fillet_edges`]                       | `SmApiFilletEdges`                                    |
//!
//! Deviations:
//! - The stages are the session's (see [`FilletSession`]); `do_fillet` runs them in
//!   `DoFilleting`'s order: `CreateFilletBrep`, the fillet-face check, `ShrinkGeometry`,
//!   then (no piecewise merge: `SmApiFilletEdges` never sets it) the edge map and topology
//!   insertion, the global merge, the trim/delete/merge and this module's clean-up.
//! - A curve the stitch copies onto the target takes the target's curve identity of its
//!   edge (`CacheObject::Curve`), as the B-rep's own stitch does.
//! - No recursion.

use super::{invalid, FilletBody, FilletObject, FilletSession};
use crate::blend::fillet_executive::p3::ErrorInfo;
use crate::error::{ignore_status, EntityRef, Result, SolidError};
use crate::geom::object::NotifyOperation;
use crate::kernel::Kernel;
use crate::topology::{Brep, EdgeId};

const OP: &str = "blend::fillet_executive::session::finish";

/// `SM_FILERR_BAD_MERGE`.
const FILERR_BAD_MERGE: i32 = 7;

impl FilletSession<'_> {
    /// `DoFilleting`'s global merge (`m_bDoGlobalMerge`): `IntersectInsertRelate` of
    /// the target and the intermediate fillet body at the target's approximation
    /// tolerance. Nothing when the global merge is off.
    ///
    /// # Errors
    /// Missing insertion state, or those of `IntersectInsertRelate` (its mutations kept).
    pub fn global_merge(&mut self, k: &mut Kernel) -> Result<()> {
        if !self.executive.options.global_merge {
            return Ok(());
        }
        // GetApproxTol3d of the target (SmTol::GetApproxTol3d(SmBrep*))
        let tolerance = crate::geom::tol::approx_tol3d(k, 0.5 * self.objects.target.tolerance());
        self.executive.intersector.set_this_approx_tol_3d(tolerance);
        self.with_insertion_geometry(k, |k, intersector, geometry| {
            intersector
                .intersect_insert_relate(k, geometry, FilletObject::Intersector)
                .map(|_| ())
        })
    }

    /// `DoFilleting` after the fillet faces are merged into the target: clears the
    /// intersector's fillet body and deletes it, stitches the target (`SmStitch` with
    /// region nesting and `m_bIgnoreProblems`, up to 5 tries doubling the tolerance from
    /// the target's), removes the topological edges the merge tagged (the first user of
    /// each retained merge attribute), runs `BooleanPostProcess(TRUE)`, notifies every
    /// target face `SM_NO_PRE_EDIT` and orients the trimmed surfaces (statuses ignored,
    /// as upstream). The merge attributes are destroyed on every return.
    ///
    /// # Errors
    /// [`SolidError::InvalidInput`] before the merge ran; `SM_ERR` (with the executive's
    /// error noted as a bad merge) when the fifth stitch still fails; those of
    /// `RemoveTopoEdges` and `BooleanPostProcess`.
    pub fn finish_fillet_merge(&mut self, k: &mut Kernel) -> Result<()> {
        let attributes = self
            .objects
            .merge_attributes
            .take()
            .ok_or_else(|| invalid("fillet faces have not been merged"))?;
        let result = self.finish_with_attributes(k, &attributes);
        let mut cleanup = Ok(());
        for attribute in attributes {
            let destroyed = self.objects.attributes.destroy(attribute);
            if cleanup.is_ok() {
                cleanup = destroyed;
            }
        }
        result.and(cleanup)
    }

    fn finish_with_attributes(
        &mut self,
        k: &mut Kernel,
        attributes: &[crate::geom::attribute::AttributeHandle],
    ) -> Result<()> {
        let fillet = Some(FilletObject::IntermediateBrep);
        if self.executive.intersector.other_brep() == fillet {
            self.executive.intersector.set_other_brep(None);
        }
        if self.executive.intersector.primary_brep() == fillet {
            self.executive.intersector.set_primary_brep(None);
        }
        // delete m_pFilletBrep
        self.objects.intermediate = None;
        self.executive.fillet = None;

        let editing = self.objects.target.editing_enabled();
        self.objects.target.set_editing_enabled(true);
        let result = self.clean_up_target(k, attributes);
        self.objects.target.set_editing_enabled(editing);
        result
    }

    fn clean_up_target(
        &mut self,
        k: &mut Kernel,
        attributes: &[crate::geom::attribute::AttributeHandle],
    ) -> Result<()> {
        // stitch the merged fillet faces to the original faces
        let mut tolerance = self.objects.target.tolerance();
        for attempt in 0..5 {
            let lamina = self.stitch_target(k, tolerance);
            match lamina {
                Ok(0) => break,
                Ok(_) => {}
                Err(error) if error.is_unported() => return Err(error),
                Err(_) if attempt == 4 => {
                    self.executive.error.code = FILERR_BAD_MERGE;
                    self.executive.error.comment =
                        "Fillet Stitching Failure: Error Stitching Fillets To Original Brep".into();
                    return Err(SolidError::sm_err(OP));
                }
                Err(_) => {}
            }
            tolerance *= 2.0;
        }
        self.executive
            .intersector
            .set_primary_brep(Some(FilletObject::Body(FilletBody::Target)));

        // the topological edges the merge tagged: the first user of each attribute
        let serial = self.objects.target_serial;
        let mut topo_edges: Vec<Option<EdgeId>> = Vec::new();
        for &attribute in attributes {
            let Some(&user) = self.objects.attributes.users(attribute)?.first() else {
                continue;
            };
            let edge = self.objects.target.edges().iter().map(|(e, _)| e).find(|&e| {
                crate::topology::brep::p3::topology_ref(serial, EntityRef::Edge(e)).key == user
            });
            topo_edges.push(edge);
        }
        crate::blend::fillet_executive::p1::remove_topo_edges(k, self.objects.target, &topo_edges)?;
        crate::boolean::merge::p1::boolean_post_process(
            k,
            self.objects.target,
            &crate::boolean::merge::p1::CleanupRelations::default(),
            true,
        )?;

        // the faces' Notify(SM_NO_PRE_EDIT) and OrientTrimmedSurfaces (statuses ignored)
        let target = &mut *self.objects.target;
        let key_of = |entity| crate::topology::brep::p3::topology_ref(serial, entity);
        crate::topology::brep::p3::seed_records(
            &mut self.objects.notification_records,
            target,
            &key_of,
        );
        let faces = crate::topology::brep::p1::faces(target, None)?;
        // (as the B-rep's own StitchAndOrient: the orientation notifies through the
        // caller's environment, with the kernel context kept)
        let mut env = crate::topology::edge::NotifyEnv {
            records: &mut self.objects.notification_records,
            attributes: &mut self.objects.attributes,
            context: None,
            key_of: &key_of,
        };
        for &face in &faces {
            let caller = key_of(EntityRef::Face(face));
            ignore_status(crate::topology::face::p3::notify(
                target,
                face,
                &mut env,
                None,
                NotifyOperation::PreEdit,
                [Some(caller), None, None],
            ))?;
        }
        ignore_status(crate::topology::brep::p3::orient_trimmed_surfaces(
            k,
            target,
            &mut env,
            &mut crate::topology::topology_traverser::Marks::default(),
            true,
            false,
            None,
        ))
    }

    /// One `SmStitch::DoStitching` pass over the target (region nesting, ignoring
    /// glue problems); the number of lamina edges left.
    fn stitch_target(&mut self, k: &mut Kernel, tolerance: f64) -> Result<u64> {
        let serial = self.objects.target_serial;
        let target = &mut *self.objects.target;
        let key_of = |entity| crate::topology::brep::p3::topology_ref(serial, entity);
        let keys = |object| match object {
            crate::topology::brep_cache::CacheObject::Brep => {
                key_of(EntityRef::Brep(crate::topology::BrepId::from_index(0)))
            }
            crate::topology::brep_cache::CacheObject::Other(owner) => owner,
            object => crate::topology::brep::p2::caller_cache_key(&key_of, object),
        };
        let owner_of = |_| None;
        let identity = crate::query::brep::BrepIdentity {
            brep_id: crate::topology::BrepId::from_index(
                u32::try_from(serial)
                    .map_err(|_| invalid("target body identity exceeds query domain"))?,
            ),
            keys: &keys,
            owner_of: &owner_of,
        };
        crate::topology::brep::p3::seed_records(
            &mut self.objects.notification_records,
            target,
            &key_of,
        );
        crate::topology::brep::p3::with_lent_context(k, |k, context| {
            let mut env = crate::topology::edge::NotifyEnv {
                records: &mut self.objects.notification_records,
                attributes: &mut self.objects.attributes,
                context: Some(context),
                key_of: &key_of,
            };
            // a copied stitch curve takes its edge's target curve identity
            let curve_ref = |edge: EdgeId, _: i64| {
                Ok(keys(crate::topology::brep_cache::CacheObject::Curve(
                    edge,
                )))
            };
            let mut edits = crate::heal::stitch::p1::BorrowedStitchEdits {
                env: &mut env,
                curve_ref: &curve_ref,
            };
            let mut callback = crate::heal::stitch::p1::DefaultStitchCallback::new();
            let mut stitch =
                crate::heal::stitch::p1::Stitch::new(&mut callback, tolerance, true, true);
            stitch.ignore_problems = true;
            stitch.query_identity = Some(identity);
            stitch.edits = Some(&mut edits as &mut dyn crate::heal::stitch::p1::StitchEditContext);
            stitch
                .do_stitching(
                    k,
                    target,
                    None,
                    None,
                    crate::heal::stitch::p2::MinGapRequest::default(),
                )
                .map(|stitched| stitched.lamina_edges)
        })?
    }

    /// `SmFilletExecutive::DoFillet` after `SetFilletParameters`: `CreateFilletCorners`,
    /// then [`FilletSession::do_filleting`].
    ///
    /// # Errors
    /// The first failing stage's error; the mutations of the stages before it are kept.
    pub fn do_fillet(
        &mut self,
        k: &mut Kernel,
        mut error_info: Option<&mut ErrorInfo<FilletObject>>,
    ) -> Result<()> {
        self.create_corners(k, error_info.as_deref_mut())?;
        self.do_filleting(k, error_info)
    }

    /// `SmFilletExecutive::DoFilleting` (under the filleting context flags):
    /// `CreateFilletBrep` (a failure or an empty fillet body is `SM_ERR`),
    /// `ShrinkGeometry` of the fillet faces, the edge map and topology insertion, the
    /// global merge, the trim, delete and merge of the original faces and
    /// [`FilletSession::finish_fillet_merge`].
    ///
    /// # Errors
    /// The first failing stage's error; the mutations of the stages before it are kept.
    pub fn do_filleting(
        &mut self,
        k: &mut Kernel,
        error_info: Option<&mut ErrorInfo<FilletObject>>,
    ) -> Result<()> {
        self.with_filleting_context(k, |session, k| {
            // (FILEXEC_ERR: any failure is SM_ERR; a missing port still defers)
            match session.create_fillet_brep(k, error_info) {
                Ok(()) => {}
                Err(error) if error.is_unported() => return Err(error),
                Err(_) => return Err(SolidError::sm_err(OP)),
            }
            let empty = session
                .intermediate_brep()
                .map_or(Ok(true), |body| {
                    crate::topology::brep::p1::faces(body, None).map(|faces| faces.is_empty())
                })?;
            if empty {
                return Err(SolidError::sm_err(OP));
            }
            session.shrink_fillet_faces(k)?;
            session.insert_fillet_intersection_topology(k)?;
            session.global_merge(k)?;
            let mut deleted = super::original_trim::OriginalSurfaceDeletion::default();
            session.trim_delete_and_merge_fillet_faces(k, &mut deleted)?;
            session.finish_fillet_merge(k)
        })
    }
}

/// `SmApiFilletEdges(pBrep, rEdges, dRadius, lXSectType, lContinuity, dThumbweight)`:
/// constant-radius fillets of `edges` of `brep` (cross section `0` linear, `2` blend
/// curve, otherwise circular), with the global merge, then `StitchAndOrient`.
///
/// # Errors
/// `SM_ERR_INVALID_INPUT` for no edges; those of the fillet stages and of the final
/// stitch and orientation.
pub fn fillet_edges(
    k: &mut Kernel,
    brep: &mut Brep,
    edges: &[EdgeId],
    radius: f64,
    xsect_type: u32,
    continuity: u32,
    thumbweight: f64,
) -> Result<()> {
    if edges.is_empty() {
        return Err(SolidError::Upstream {
            operation: OP,
            code: crate::error::SM_ERR_INVALID_INPUT,
            entity: None,
        });
    }
    {
        let mut session = FilletSession::new(k, brep)?;
        session.prepare_constant_radius_edges(
            k,
            edges,
            radius,
            xsect_type,
            continuity,
            thumbweight,
        )?;
        session.do_fillet(k, None)?;
    }
    // pBrep->StitchAndOrient() (status ignored upstream)
    ignore_status(crate::topology::brep::p5::stitch_and_orient(
        k,
        brep,
        &crate::topology::brep::p5::StitchAndOrientOptions::default(),
    ))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::math::Point3;

    fn boxed(k: &mut Kernel) -> Brep {
        crate::primitive::api_primitives::create_box(
            k,
            Point3::new(0.0, 0.0, 0.0),
            10.0,
            20.0,
            30.0,
        )
        .unwrap()
    }

    #[test]
    fn fillet_edges_needs_an_edge() {
        let mut k = Kernel::default();
        let mut b = boxed(&mut k);
        let r = fillet_edges(&mut k, &mut b, &[], 1.0, 1, 1, 1.0);
        assert!(matches!(
            r,
            Err(SolidError::Upstream { code, .. }) if code == crate::error::SM_ERR_INVALID_INPUT
        ));
    }

    /// One box edge rounded end to end: the six faces plus the fillet face, 15 edges
    /// and 10 vertices (the single-edge box counts the C++ reference reports).
    #[test]
    fn fillet_one_box_edge_end_to_end() {
        let mut k = Kernel::default();
        let mut b = boxed(&mut k);
        let e = b.edges().ids().next().unwrap();
        fillet_edges(&mut k, &mut b, &[e], 1.0, 1, 1, 1.0).unwrap();
        let faces = crate::topology::brep::p1::faces(&b, None).unwrap().len();
        let edges = crate::topology::brep::p1::edges(&b, None).unwrap().len();
        let vertices = crate::topology::brep::p1::vertices(&b, None).unwrap().len();
        assert_eq!((faces, edges, vertices), (7, 15, 10));
    }
}
