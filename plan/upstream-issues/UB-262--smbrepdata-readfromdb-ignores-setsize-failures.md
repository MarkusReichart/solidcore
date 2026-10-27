#### UB-262. `SmBrepData::ReadFromDB` sizes its arrays from unchecked file counts [PSIRT]
- **Severity:** Crash / undefined behaviour
- **Area:** File input
- **Location:** `source/SMLib/src/SmBrepData.cpp`, `SmBrepData::ReadFromDB` (def. `:2933`), at these array-sizing calls and fill loops:

  | array | `SetSize` | fill loop |
  |---|---|---|
  | `m_vRegions` | `:3049` | `:3052` |
  | `m_vShells` | `:3075` | `:3078` |
  | `m_vFaceuses` | `:3157` | `:3160` |
  | `m_vFaces` | `:3183` | `:3186` |
  | `m_vLoops` | `:3208` | `:3211` |
  | `m_vEdgeuses` | `:3292` | `:3295` |
  | `m_vEdges` | `:3316` | `:3319` |
  | `m_vVertices` | `:3340` | `:3343` |
  | `m_vUVCurves` | `:3363` | `:3366` |
  | `m_v3DCurves` | `:3400` | `:3403` |
  | `m_vSurfaces` | `:3437` | `:3438` |
- **What happens:** `ReadFromDB` sizes its arrays from counts read from the file without checking `SetSize`'s result, then fills them. In release builds the array index is unchecked, so a crafted file can cause out-of-bounds writes.
- **How found:** tangency-labs/solidcore#845 (source reading of 1.0.1). Re-read here at the pin.
- **Confirmed:** no (read at c9c5979 `source/SMLib/src/SmBrepData.cpp:3049`-`3438`).
- **Rust port:** not affected. `BrepData::read_from_db` (`crates/solidcore/src/topology/brep_data.rs`) checks every count against the remaining input before reading any record (`checked_count`, `:309`-`316`). It then appends records one at a time.
- tangency-labs/solidcore#845.
