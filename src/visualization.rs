use minifb::{Key, Window, WindowOptions};
use crate::{BIOME_COLORS, Game, MAP_SIZE};

pub fn visualize(game: &Game) {
    let mut window = Window::new(
        "Biome Map Generator",
        MAP_SIZE,
        MAP_SIZE,
        WindowOptions::default(),
    ).expect("Nie udało się otworzyć okna");

    // Bufor dla okna (musi być u32 - format 0xRRGGBB)
    let mut buffer: Vec<u32> = vec![0; MAP_SIZE * MAP_SIZE];

    while window.is_open() && !window.is_key_down(Key::Escape) {
        // Konwersja Biome -> Kolor
        // Możesz to też zrównoleglić przez Rayon!
        for (i, row) in game.board.iter().enumerate() {
            for (j, &cell) in row.iter().enumerate() {
                buffer[i * MAP_SIZE + j] = BIOME_COLORS[cell as usize];
            }
        }

        window.update_with_buffer(&buffer, MAP_SIZE, MAP_SIZE).unwrap();
    }
}