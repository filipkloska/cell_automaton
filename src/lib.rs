pub mod cell;
pub mod game;
pub mod visualization;

pub use crate::game::Game;
pub use crate::cell::Biome;

pub use::colored::Colorize;
pub use minifb::{Window, WindowOptions, Key};

pub const READ_WINDOW_WIDTH: usize = 2;
pub const READ_WINDOW_HEIGHT: usize = 2;
pub const SAFETY_OFFSET_SIZE: usize = 2;
pub const MAP_SIZE: usize = 1024;
pub const OCEAN_BIOME_CHANCE: [f64; 6] = [0.8, 0.2, 0.1, 0.0, -0.1, -0.1]; //ocean, sea, coast, lowlands
pub const SEA_BIOME_CHANCE: [f64; 6] = [0.2, 0.2, 0.8, 0.0, -0.1, -0.1];
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


const BIOME_COLORS: [u32; 8] = [
    0x000000, // Empty (Czarny)
    0x0000FF, // Ocean (Niebieski)
    0x4444FF, // Sea (Jasnoniebieski)
    0xFFFFCC, // Coast (Piaskowy)
    0x228B22, // Lowlands (Zielony)
    0x006400, // Highlands (Ciemnozielony)
    0x808080, // Mountains (Szary)
    0xFF00FF, // DUMMY (Magenta - dla testu granic)
];



