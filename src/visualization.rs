use minifb::{Key, Window, WindowOptions};
use crate::{BIOME_COLORS, Game, MAP_SIZE};

pub fn visualize<const N: usize>(game: &Game<N>) {
    let mut window = Window::new(
        "Biome Map Generator",
        game.map_size,
        game.map_size,
        WindowOptions::default(),
    ).expect("Nie udało się otworzyć okna");

    // Bufor dla okna (musi być u32 - format 0xRRGGBB)
    let mut buffer: Vec<u32> = vec![0; MAP_SIZE * MAP_SIZE];

    while window.is_open() && !window.is_key_down(Key::Escape) {
        for (i, row) in game.board.iter().enumerate() {
            for (j, &cell) in row.iter().enumerate() {
                buffer[i * MAP_SIZE + j] = BIOME_COLORS[cell as usize];
            }
        }

        window.update_with_buffer(&buffer, MAP_SIZE, MAP_SIZE).unwrap();
    }
}