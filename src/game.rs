use std::{collections::VecDeque, usize};

use crate::{Cell, Biome, State, MAP_SIZE, Colorize};
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
        
        
        let start_idx = MAP_SIZE/2;
        self.board[start_idx][start_idx].transform_cell(Biome::Ocean);
        let mut transform_indices = 
            VecDeque::from([(start_idx,start_idx)]);
        
        /*add to queue cells to transform
            find neighbours of each 
            neighbours = (pop_front.find)
            next calculate values from all 8(calculate)
            transform i-1 j, i+1 j, i j-1, i j+1
            add these to queue
            repeat 
            */
        // self.print_board_biome();
        // self.print_board_state();
        // while let Some(current) = transform_indices.pop_front() {
        //     println!("{:?}", transform_indices);
        //     let neighbours = self.find_neighbours(current);
        //     let biome = self.calculate_neighbours(neighbours);
        //     self.board[current.0][current.1].transform_cell(biome);
            
        //     let newly = self.find_adjacent_cells(current);
        //     for n in newly {
        //         transform_indices.push_back(n);
        //     }
        //     self.print_board_biome();
        //     self.print_board_state();   
        // }
        while !transform_indices.is_empty() {
            let mut next_layer = Vec::new();
            let layer_len = transform_indices.len();
            for _ in 0..layer_len {
                let current = match transform_indices.pop_front() {
                    Some(val) => val,
                    None => panic!("Popped empty queue")
                };
                let neighbours = self.find_neighbours(current);
                let b = self.calculate_neighbours(neighbours);
                
                self.board[current.0][current.1].transform_cell(b);
                
                let adjencent = self.find_adjacent_cells(current);
                for indices in adjencent {
                    next_layer.push(indices);
                }
            }
            for &(i, j) in &next_layer {
                self.board[i][j].state = State::Buffered;
            }

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
            if let (Some(ni), Some(nj)) = (i.checked_add_signed(oi), j.checked_add_signed(oj)) {
                if ni < MAP_SIZE && nj < MAP_SIZE {
                    match self.board[ni][nj].state {
                        State::Empty => {
                            self.board[ni][nj].state = State::Buffered;
                            adjacent.push((ni, nj));
                        }
                    _ => continue, // Buffered lub Transformed → pomijamy
                    }
                }
            }
        }

    adjacent
}


    /*Takes a queue of cell indices to calculate neighbours of each cell and transforms the board */
    /*for now a place holder */
    fn calculate_neighbours(&self, neighbours: Vec<(usize,usize)>) -> Biome {
        return Biome::Ocean;
    }

    /*Looks for neighbours for the cell*/
    fn find_neighbours(&self, (i, j): (usize, usize)) -> Vec<(usize, usize)> {
    let mut neighbours = Vec::new();
        for di in -1..=1 {
            for dj in -1..=1 {
                if di == 0 && dj == 0 {
                    continue;
                }
                if let (Some(ni), Some(nj)) = 
                    (i.checked_add_signed(di), 
                    j.checked_add_signed(dj)) {
                    if ni < MAP_SIZE && nj < MAP_SIZE {
                        neighbours.push((ni, nj));
                    }
                }
            }
        }

    neighbours
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
                    Biome::Sea => print!("{}"," S ".bright_blue()),
                    Biome::Coast => print!("{}"," C ".bright_yellow()),
                    Biome::Highlands => println!("{}"," H ".yellow()),
                    Biome::Lowlands => println!("{}"," L ".bright_green()),
                    Biome::Mountains => println!("{}"," M ".green())
                }
            }
            println!()
        }
    }
}