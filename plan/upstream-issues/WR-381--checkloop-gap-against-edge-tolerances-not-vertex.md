#### WR-381. `SmTrimmingTools::CheckLoop` measures the gap at a vertex against the edges' tolerances, not the vertex's
- **Severity:** Wrong result (a check failure on a body whose tolerances are consistent)
- **Area:** Trimming / checking
- **Location:** `source/SMLib/src/SmTrimmingTools.cpp` (around `:3127`), `SmTrimmingTools::CheckLoop`: `SmXSectTol3d sXSectTol3d = SmTol::GetXSectTol3d(sEdgeuses[ii]->GetEdge(), sEdgeuses[previous]->GetEdge())`, then `dEndStartDist <= sXSectTol3d` (the source comment notes it is tighter than the tolerance model).
- **What happens:** where two edges meet at a vertex whose tolerance was widened to cover the gap, the check still compares the gap with the two edges' tolerances and fails (closure 3, a `SM_LIST_3` report).
- **Trigger:** a union of two cylinders where an intersection curve ends 1.3e-5 short of its vertex: edge tolerances 5.8e-6 (intersection tolerance 1.16e-5), vertex tolerance 2.6e-5 (tangency-labs/solidcore#476, SolidCraft regression rig on SMLib).
- **How found:** tangency-labs/solidcore#476.
- **Confirmed:** source; the port reproduces upstream with the default.
- **Fix:** the gap is within tolerance when it is within the larger of the edges' intersection tolerance and the zone tolerance of the vertex where the two edge uses meet.
- **Rust port:** opt-in `Kernel::fixes.check_loop_vertex_tolerance` (default `false`, upstream's test), read by `loop_gap_tolerance` (`check_loop`'s gap tolerance) in `crates/solidcore/src/topology/trimming_tools.rs`; documented in its Deviations. Test `issue_476_check_loop_vertex_tolerance`. With the default, `topology.trimming_tools` (which runs `check_loop`) compares as before. Fixes tangency-labs/solidcore#476.
