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

