use crate::{Cell, Biome, State, MAP_SIZE, Colorize};
pub struct Game {
    board: [[Cell; MAP_SIZE]; MAP_SIZE]
}

impl Game {
    pub fn new() -> Game {
        // initialize board
        let board = 
        [[Cell{ biome: Biome::Empty, state: State::Empty}; MAP_SIZE]; MAP_SIZE];
        return Game {board: board};
    }
    
    pub fn run(&mut self) {
        let mut step = 1;
        let start_idx = MAP_SIZE/2;
        // while step < MAP_SIZE * MAP_SIZE {
            


        // }
        self.board[start_idx][start_idx].transform_cell(Biome::Ocean);
        self.print_board_biome();
    }

    fn calculate_neighbours(&self) {

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