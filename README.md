# Running solidcore's difftest against the SMLib fork

Local scripts, not a finished tool. They run `crates/difftest` with the C++ side built from
tangency-labs/smlib instead of the pinned v1.0.1. Nothing here is claimed as work in progress:
take whatever is useful. Last run: solidcore b0e925d6c against the fork at 5e8c9c1, Linux.

## What `run.sh` expects

Paths are hard-coded for one machine; change `B`, `E`, `R`, `W` and `F` at the top.

- A solidcore checkout in `$B/repo` that holds the pinned `upstream/` (`TestFiles`, `tests`).
- The fork in `$F`, and its `libSMLib.a` and `libSM_API.a` built with the harness's CMake project
  (`-DSMLIB_UPSTREAM=$F -DSC_BUILD_TESTS=OFF`), linked from `$SOLIDCORE_CPP_DIR/lib`, next to a
  hand-written `ffi.env` with `INCLUDE_DIRS`, `DEFINES`, `CXXFLAGS`, `LIBS` and `SYSLIBS`.

## What it does

1. A worktree of solidcore `main`, with `upstream/source` pointing at the fork.
2. A patch to `Kernel::new`: every flag on, `FORKDIFF_NONE=1` starts from none, `FORKDIFF_ON` and
   `FORKDIFF_OFF` name flags by name.
3. `sync_copies.py`: file-static functions the shims copied from v1.0.1 become the fork's, where
   the copy is pure. `stale_copies.py` lists the rest.
4. Shims that don't compile against the fork leave. Their entry points become aborting stubs
   (`stubs.py`, `missing.txt`), and the build links with `--unresolved-symbols=ignore-all`, so
   those units crash and report nothing. `restore_small.py` adapts five shims by hand.
5. `difftest --green` with all flags on, then off. `compare.py` sorts every case into a class.

`isolate.py` finds the flags behind a case by delta debugging. `known.py` takes out the two known
false alarms. `triage.py` sorts failures by kind.

## Known gaps

- `curve_composite_curve_p1.cpp` copies `sm_CreateCorner` with the `ComputeCrvsNormal` status
  ignored, where the fork has `SER`. `run.sh` doesn't patch that, which gives 33 false refusals.
- `surface_cone_p2.cpp` and `surface_surface_p1.cpp` leave only because their copies call changed
  signatures. They could be brought along.
- The units of the fillet solver, geometry and corners don't report at all.
