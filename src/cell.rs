#[repr(u8)]
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Biome{
    Empty = 0,
    Ocean = 1,
    Sea = 2,
    Coast = 3,
    Lowlands = 4,
    Highlands = 5,
    Mountains = 6,
    Checked = 7
}

impl Biome {
    pub fn from_u8(value: u8) -> Biome {
        match value {
            0 => Biome::Empty,
            1 => Biome::Ocean,
            2 => Biome::Sea,
            3 => Biome::Coast,
            4 => Biome::Lowlands,
            5 => Biome::Highlands,
            6 => Biome::Mountains,
            7 => Biome::Checked,
            _ => Biome::Empty,
        }
    }
}

