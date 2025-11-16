pub mod biome;
pub mod cell;
pub mod game;

pub const MAP_SIZE: usize = 20;
pub const OCEAN_BIOME_CHANCE: [f64; 6] = [0.7, 0.2, 0.1, 0.0, 0.0, 0.0]; //ocean, sea, coast
pub const SEA_BIOME_CHANCE: [f64; 6] = [0.2, 0.2, 0.6, 0.0, 0.0, 0.0]; //ocean, sea, coast
pub const COAST_BIOME_CHANCE: [f64; 6] = [0.1, 0.2, 0.7, 0.0, 0.0, 0.0]; //ocean, sea, coast, lowlands


pub use crate::{cell::Cell, biome::Biome, cell::State, game::Game};
pub use::colored::Colorize;