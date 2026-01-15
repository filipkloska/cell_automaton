use cell_automaton::{Game, visualization};
fn main() {
    let mut game = Game::new();
    game.run();
    visualization::visualize(&game);
}
