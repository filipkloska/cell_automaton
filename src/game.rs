use std::{collections::VecDeque, usize};
use rand::Rng;
use crate::{Biome, Cell, Colorize, MAP_SIZE, 
    OCEAN_BIOME_CHANCE, SEA_BIOME_CHANCE, COAST_BIOME_CHANCE, 
    State};
pub struct Game {
    board: [[Cell; MAP_SIZE]; MAP_SIZE]
}

impl Game {
    /*Constructor initializing the board */
    pub fn new() -> Game {
        let board = 
        [[Cell{ biome: Biome::Empty, state: State::Empty}; MAP_SIZE]; MAP_SIZE];
        return Game {board: board};
    }
    /*Basic game loop */
    
    pub fn run(&mut self) {
        let mut densities: [f64; 6] = [0.99; 6];
        let start_idx = MAP_SIZE/2;
        self.board[start_idx][start_idx].transform_cell(Biome::Ocean);
        let mut transform_indices = 
            VecDeque::from([(start_idx,start_idx)]);
        
        
        while !transform_indices.is_empty() {
            //next layer takes adjencent cells
            let mut next_layer = Vec::new();
            let layer_len = transform_indices.len();
            for _ in 0..layer_len {
                let current = match transform_indices.pop_front() {
                    Some(val) => val,
                    None => panic!("Popped empty queue")
                };
                let adjencent = self.find_adjacent_cells(current);
                for indices in &adjencent {
                    let around = self.find_cells_around((indices.0, indices.1));
                    
                    self.board[indices.0][indices.1].
                    transform_cell(self.calculate_biome_for_cell(current,around, &mut densities));
                }
                
                for indices in adjencent {
                    next_layer.push(indices);
                }
            }
            //i think this will work for multithreading
            // for &(i, j) in &next_layer {
            //     self.board[i][j].state = State::Buffered;
            // }

            for n in next_layer {
                transform_indices.push_back(n);
            }
            self.print_board_biome();
        }
        
    }
    


    fn find_adjacent_cells(&mut self, (i, j): (usize, usize)) -> Vec<(usize, usize)> {
        let mut adjacent = Vec::new();
        let offsets = [(-1, 0), (1, 0), (0, -1), (0, 1)];

        for (oi, oj) in offsets {
            if let (Some(ni), Some(nj)) = 
            (i.checked_add_signed(oi), j.checked_add_signed(oj)) {
                if ni < MAP_SIZE && nj < MAP_SIZE {
                    if let State::Empty = self.board[ni][nj].state {
                        self.board[ni][nj].state = State::Buffered;
                        adjacent.push((ni, nj));
                    }
                }
            }
        }

    return adjacent;
}


    /*Takes a queue of cell indices to calculate neighbours of each cell and transforms the board */
    /*for now a place holder */
    fn calculate_biome_for_cell(&self, cell: (usize, usize), neighbours: Vec<(usize,usize)>, densities: &mut [f64; 6]) -> Biome {
        let mut chances: [f64; 6] = [0.0; 6];
        
        let mut mult = 1.0;

        for (i, j) in neighbours {
            
            if self.board[i][j].biome == Biome::Empty {
                continue;
            }
            if i != cell.0 && j != cell.1 {
                mult = 0.5;
            }
            else {
                mult = 1.0;
            }
            
            match self.board[i][j].biome {
                Biome::Ocean => {
                    for i in 0..6 {
                        chances[i] += OCEAN_BIOME_CHANCE[i] * mult * densities[i];
                    }
                    //densities[0] *= 0.99;
                },
                Biome::Sea => {
                    for i in 0..6 {
                        chances[i] += SEA_BIOME_CHANCE[i] * mult * densities[i];
                    }
                    //densities[1] *= 0.99;
                },
                Biome::Coast => {
                    for i in 0..6 {
                        chances[i] += COAST_BIOME_CHANCE[i] * mult * densities[i];
                    }
                    //densities[2] *= 0.99;
                },
                _ => {}
            }
        } 
        let mut sum = 0.0;
        
        for val in &chances {
            sum += *val;
        }
        for val in &mut chances {
            *val /= sum;
        }
        
        let probability: f64 = rand::rng().random_range(0..101) as f64 / 100.0;
        if probability <= chances[0] {
            Biome::Ocean
        }
        else if probability <= chances[0] + chances[1] {
            Biome::Sea
        }
        else {
            Biome::Coast
        }
    }

    /*Looks for neighbours around the cell*/
    fn find_cells_around(&self, (i, j): (usize, usize)) -> Vec<(usize, usize)> {
    let mut neighbours = Vec::new();
        for di in -1..=1 {
            for dj in -1..=1 {
                if di == 0 && dj == 0 {
                    continue;
                }
                if let (Some(ni), Some(nj)) = 
                (i.checked_add_signed(di), j.checked_add_signed(dj)) {
                    if ni < MAP_SIZE && nj < MAP_SIZE {
                        neighbours.push((ni, nj));
                    }
                }
            }
        }

    return neighbours;
}
    pub fn print_board_state(&self) {
        for row in self.board{
            for cell in row {
                match cell.state {
                    State::Empty => print!("{}"," . "),
                    State::Transformed => print!("{}"," T "),
                    State::Buffered => print!("{}"," B "),
                }
            }
            println!();
        }
    }

    pub fn print_board_biome(&self) {
        for row in self.board{
            for cell in row {
                match cell.biome {
                    Biome::Empty => print!("{}"," . ".bright_black()),
                    Biome::Ocean => print!("{}"," O ".blue()),
                    Biome::Sea => print!("{}"," S ".magenta()),
                    Biome::Coast => print!("{}"," C ".bright_yellow()),
                    Biome::Highlands => println!("{}"," H ".yellow()),
                    Biome::Lowlands => println!("{}"," L ".bright_green()),
                    Biome::Mountains => println!("{}"," M ".green())
                }
            }
            println!()
        }
        println!()
    }
}