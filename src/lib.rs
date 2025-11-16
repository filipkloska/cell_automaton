pub mod biome;
pub mod cell;
pub mod game;

pub const MAP_SIZE: usize = 20;
pub const OCEAN_BIOME_CHANCE: [f64; 6] = [0.4, 0.2, 0.1, 0.0, -0.1, -0.1]; //ocean, sea, coast, lowlands
pub const SEA_BIOME_CHANCE: [f64; 6] = [0.2, 0.2, 0.6, 0.0, -0.1, -0.1];
pub const COAST_BIOME_CHANCE: [f64; 6] = [0.1, 0.2, 0.1, 0.6, 0.0, 0.0];
pub const LOWLANDS_BIOME_CHANCE: [f64; 6] = [0.0, 0.0, 0.1, 0.5, 0.2, 0.1];
pub const HIGHLANDS_BIOME_CHANCE: [f64; 6] = [-0.1, -0.1, 0.0, 0.4, 0.4, 0.2];
pub const MOUNTAINS_BIOME_CHANCE: [f64; 6] = [-0.1, -0.1, 0.0, 0.1, 0.7, 0.2];

// pub const OCEAN_BIOME_CHANCE: [f64; 6] = [ 1.0,  0.7,  0.4, -0.2, -0.8, -1.0];
// pub const SEA_BIOME_CHANCE:   [f64; 6] = [ 0.7,  1.0,  0.8, -0.1, -0.6, -0.9];
// pub const COAST_BIOME_CHANCE: [f64; 6] = [ 0.4,  0.8,  1.0,  0.5, -0.2, -0.5];
// pub const LOWLANDS_BIOME_CHANCE:[f64;6] =[-0.2, -0.1,  0.5,  1.0,  0.6,  0.2];
// pub const HIGHLANDS_BIOME_CHANCE:[f64;6] =[-0.8, -0.6, -0.2,  0.6,  1.0,  0.8];
// pub const MOUNTAINS_BIOME_CHANCE:[f64;6] =[-1.0, -0.9, -0.5,  0.2,  0.8,  1.0];


pub use crate::{cell::Cell, biome::Biome, cell::State, game::Game};
pub use::colored::Colorize;