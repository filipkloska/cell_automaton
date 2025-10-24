use crate::{Biome};

#[derive(Copy, Clone, Debug)]
pub struct Cell {
    pub biome: Biome,
    pub state: State   
}

//not sure if we need this, we can store fresh cells in an array in game
#[derive(Copy, Clone, Debug)]
pub enum State{
    Empty,
    Fresh,
    Old,
}

impl Cell {
    pub fn transform_cell(&mut self, biome: Biome) {
        self.biome = biome;
    }
}
