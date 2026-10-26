# topology.brep_p3: the two Mouse TurnToNURBS cases notify the trim and stretch (2026-10-25)

Updates: plan/findings/topology.brep_p3--turntonurbs-uv-trim-curve-edits-notify-2026-10-25.md

- The stretch was already ported, behind `Context::trim_uv_curves_at_corner_drops`
  (pinned-`c9c5979` toggle, default `false` for WR-255 / #570), but without notifications.
- Now `face::p2::fix_uv_trim_curves_notified` sends `SM_NO_PRE_EDIT` /
  `SM_NO_POST_EDIT` around every `Trim` and `EditParameterization` of the stretch that
  changes the curve's interval (`notified_interval_edit`). Upstream returns before its
  `PRE_EDIT` exactly when the interval would not change. Caller identity and owner
  links are as for `EditEndPoint`.
- The brep_p3 runner sets the toggle, as topology.brep_p10's runner does.
- `crates/difftest/known_mismatches/topology.brep_p3.toml` is deleted.
- Not decided here: under #1018 (defaults keep SMLib), the default should arguably be
  the stretch, with WR-255 as a `k.fixes` opt-in.
