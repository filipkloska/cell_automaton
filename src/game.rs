use std::{collections::VecDeque};

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
    pub fn run(&mut self, ) {
        let start_idx = MAP_SIZE/2;
        let mut transform_indices = VecDeque::from([(start_idx,start_idx)]);
        
        self.board[start_idx][start_idx].transform_cell(Biome::Ocean);
        for _ in 1..=MAP_SIZE*MAP_SIZE {
            
        }
        self.print_board_biome();
    }


    /*Takes a queue of cell indices to calculate neighbours of each cell and transforms the board */
    fn calculate_neighbours(&self, queue: &mut VecDeque<(usize,usize)>) {
        for cell in queue {
            
        }
    }

    /*Looks for neighbours for the cell*/
    fn find_neighbours(&self, cell_idx: (usize,usize)) -> Vec<(usize,usize)> {
        let neighbours = Vec::<(usize,usize)>::new();
        
        return neighbours;
    }

    pub fn print_board_state(&self) {
        for row in self.board{
            for cell in row {
                match cell.state {
                    State::Empty => print!("{}"," . "),
                    State::Fresh => print!("{}"," F "),
                    State::Old => print!("{}"," O "),
                }
            }
            println!()
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