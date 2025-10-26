pub mod biome;
pub mod cell;
pub mod game;
pub const MAP_SIZE: usize = 20;
pub use crate::{cell::Cell, biome::Biome, cell::State, game::Game};
pub use::colored::Colorize;