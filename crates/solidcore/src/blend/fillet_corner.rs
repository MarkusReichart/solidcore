//! Fillet corners: how the fillets of the edges meeting at a vertex are combined
//! (vertex blends, bevels, end caps).
//!
//! Upstream: `source/SMLib/inc/SmFilletCorner.h`, `source/SMLib/src/SmFilletCorner.cpp`.
//! The file is split into parts (`plan/modules.toml`, units `blend.fillet_corner_p1` ..
//! `blend.fillet_corner_p3`); each part is a submodule.

pub mod degenerate;
pub mod four_by_three;
pub mod mixed_3x2;
pub mod nx1_closed;
pub mod nx2_topology;
pub mod nxn_topology;
pub mod one_by_one;
pub mod open;
pub mod p1;
pub mod p2;
pub mod p3;
pub mod topology;
