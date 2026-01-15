use std::{collections::VecDeque, usize};
use rand::Rng;

use crate::{Biome, COAST_BIOME_CHANCE, HIGHLANDS_BIOME_CHANCE, LOWLANDS_BIOME_CHANCE, MAP_SIZE, MOUNTAINS_BIOME_CHANCE, OCEAN_BIOME_CHANCE, READ_WINDOW_HEIGHT, READ_WINDOW_WIDTH, SAFETY_OFFSET_SIZE, SEA_BIOME_CHANCE};

pub struct Game {
    pub board: Vec<[Biome; MAP_SIZE]>
}


#[derive(Clone, Copy)]
pub struct LimitOffsets
{
    pub read_up_offset: usize,
    pub read_down_offset: usize
}

impl Game {
    
    /*Constructor initializing the board */
    pub fn new() -> Game {
        let empty_row = [Biome::Empty; MAP_SIZE];
        let board = vec![empty_row; MAP_SIZE];
        return Game {board: board};
    }

    pub fn set_thread_count(&self) -> usize
    {
        let num_sys_threads = num_cpus::get();
        let num_threads = self.shrink_to_pow2(num_sys_threads);
        num_threads
    }

    pub fn run(&mut self){
        
        let thread_count = self.set_thread_count();
        let offsets = self.calculate_offsets(thread_count);
        let (normal_slices, safety_slices) = unsafe {
            self.setup_slices(thread_count)
        };
        std::thread::scope(|s|
        {
            for (slice, current_offset) in normal_slices.into_iter().zip(offsets.into_iter()) 
            {
                s.spawn(move ||
                {
                    Game::process_slice_worker(slice, current_offset);
                });
            }
        });
        std::thread::scope(|s| {
            for slice in safety_slices {
                let offsets = LimitOffsets {
                    read_up_offset: 0,
                    read_down_offset:0,
                };
                s.spawn(move || {
                    Game::process_slice_worker(slice, offsets);
                });
            }
        });

    }

    pub unsafe fn setup_slices(&mut self, thread_count: usize) -> (Vec<&mut [[Biome;MAP_SIZE]]>, Vec<&mut [[Biome;MAP_SIZE]]>)
    {
        let normal_band_sizes = self.calculate_normal_slice_lens(thread_count);
        let mut normal_bands = Vec::<&mut [[Biome;MAP_SIZE]]>::new();
        let mut safety_bands = Vec::<&mut [[Biome;MAP_SIZE]]>::new();
        let mut index: usize = 0;
        for i in 0..thread_count
        {
            let mut start = index;
            let mut end = index + normal_band_sizes[i];
            println!("Normal start: {}, end {}", start, end);
            unsafe 
            {
                let ptr = self.board.as_mut_ptr().add(start);
                let slice = std::slice::from_raw_parts_mut(ptr, end-start);
                normal_bands.push(slice);
            }
            start = end - 2*READ_WINDOW_HEIGHT;
            end = start + SAFETY_OFFSET_SIZE + 2*READ_WINDOW_HEIGHT;
            
            if i < thread_count-1
            {
                println!("Safety start: {}, end {}", start, end);
                unsafe 
                {
                    let ptr = self.board.as_mut_ptr().add(start);
                    let slice = std::slice::from_raw_parts_mut(ptr, end-start);
                    safety_bands.push(slice);
                }
            }
            index = end - 2*READ_WINDOW_HEIGHT
        }
        (normal_bands, safety_bands)
    }
                
    pub fn process_slice_worker(slice: &mut [[Biome; MAP_SIZE]], offsets: LimitOffsets)
    {
        let mut transform_indices = VecDeque::<(usize,usize)>::from([(slice.len()/2, MAP_SIZE/2)]);
        let mut densities: [f64; 6] = [0.99; 6];
        while !transform_indices.is_empty()
        {
            let mut next_layer = Vec::<(usize,usize)>::new();
            let indices_len = transform_indices.len();
            for _ in 0..indices_len
            {
                let current = transform_indices.pop_back().unwrap();
                let adjacent_cells = Game::get_adjacent_cell_indices(slice, current, offsets);
                for indices in adjacent_cells {
                    slice[indices.0][indices.1] = Game::calculate_biome(indices, slice, &mut densities);
                    next_layer.push(indices);
                }
            }
            for n in next_layer {
                transform_indices.push_back(n);
            }
        }
        
    }
    pub fn get_neighbours(slice: &mut [[Biome; MAP_SIZE]],index: (usize,usize)) -> Vec::<(usize,usize)>
    {
        let mut neighbours = Vec::<(usize,usize)>::new();
        let window_height = READ_WINDOW_HEIGHT as isize;
        let window_width = READ_WINDOW_WIDTH as isize;
        for wi in -window_height..=window_height
        {
            for wj in -window_width..=window_width
            {
                if wi == 0 && wj == 0
                {
                    continue;
                }
                if let (Some(ni), Some(nj)) = 
                (index.0.checked_add_signed(wi), index.1.checked_add_signed(wj))
                {
                    if ni < slice.len()
                    && nj < MAP_SIZE
                    {
                        neighbours.push((ni,nj));
                    }   
                }  
            }
        }
        neighbours
    }
    
    pub fn calculate_biome(index: (usize,usize), slice: &mut [[Biome; MAP_SIZE]], densities: &mut [f64; 6]) -> Biome
    {
        let neighbours = Game::get_neighbours(slice, index);
        let mut chances: [f64; 6] = [0.0; 6];
        let mut mult: f64;

        for (i, j) in neighbours {
            
            if slice[i][j] == Biome::Empty {
                continue;
            }
            if i != index.0 && j != index.1 {
                mult = 0.25;
            }
            else {
                mult = 1.0;
            }
            
            match slice[i][j] {
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
                Biome::Lowlands => {
                    for i in 0..6 {
                        chances[i] += LOWLANDS_BIOME_CHANCE[i] * mult * densities[i];
                    }
                    //densities[3] *= 0.99;
                },
                Biome::Highlands => {
                    for i in 0..6 {
                        chances[i] += HIGHLANDS_BIOME_CHANCE[i] * mult * densities[i];
                    }
                    //densities[4] *= 0.99;
                },
                Biome::Mountains => {
                    for i in 0..6 {
                        chances[i] += MOUNTAINS_BIOME_CHANCE[i] * mult * densities[i];
                    }
                    //densities[5] *= 0.99;
                },
                _ => {}
            }
        } 

        for val in &mut chances {
            if *val < 0.0 {
                *val = 0.0;
            }
        }
        
        let sum: f64 = chances.iter().sum();
        
        if sum == 0.0 { //If every biom is 0, all biomes are equally likely
            for val in &mut chances {
                *val = 1.0 / 6.0;
            }
        } else {
            for val in &mut chances {
                *val /= sum;
            }
        }
        let probability: f64 = rand::rng().random_range(0..100) as f64 / 100.0;
        if probability <= chances[0] {
            Biome::Ocean
        }
        else if probability <= chances[0] + chances[1] {
            Biome::Sea
        }
        else if probability <= chances[0] + chances[1] + chances[2] {
            Biome::Coast
        }
        else if probability <= chances[0] + chances[1] + chances[2] + chances[3] {
            Biome::Lowlands
        }
        else if probability <= chances[0] + chances[1] + chances[2] + chances[3] + chances[4] {
            Biome::Highlands
        }
        else if probability <= chances[0] + chances[1] + chances[2] + chances[3] + chances[4] + chances[5] {
            Biome::Mountains
        }
        else {
            Biome::Empty
        }
    }

    pub fn get_adjacent_cell_indices(slice: &mut [[Biome; MAP_SIZE]], index: (usize, usize), offset: LimitOffsets) -> Vec<(usize, usize)>
    {
        let mut adjacent_cell_indices = Vec::<(usize,usize)>::new();
        let offsets = Vec::from([(0,1), (1,0), (0,-1), (-1,0)]);
        for (oi, oj) in offsets
        {
            if let (Some(ni), Some(nj)) = 
            (index.0.checked_add_signed(oi), index.1.checked_add_signed(oj)) {
                if ni >= READ_WINDOW_HEIGHT - offset.read_up_offset 
                && ni < slice.len() - (READ_WINDOW_HEIGHT - offset.read_down_offset) 
                && nj < MAP_SIZE {
                    if let Biome::Empty = slice[ni][nj] {
                        slice[ni][nj] = Biome::Checked;
                        adjacent_cell_indices.push((ni, nj));
                    }
                }
            }
        }
        adjacent_cell_indices
    }
    pub fn calculate_offsets(&self, num_threads: usize) -> Vec<LimitOffsets>
    {
        let mut offsets=  Vec::<LimitOffsets>::new();
        for i in 0..num_threads
        {
            let mut r_u_offset = 0;
            let mut r_d_offset = 0;

            if i == 0
            {
                r_u_offset += READ_WINDOW_HEIGHT
            }
            if i == num_threads - 1
            {
                r_d_offset += READ_WINDOW_HEIGHT
            }
            let offset = LimitOffsets {
                read_down_offset: r_d_offset,
                read_up_offset: r_u_offset
            };
            offsets.push(offset);
        }
        offsets
    }

    pub fn calculate_normal_slice_lens(&self, num_threads: usize) -> Vec<usize>
    {
        let mut slice_lens = Vec::<usize>::new();
        let safety_thread_count = num_threads - 1;
        let normal_generation_count = MAP_SIZE - (safety_thread_count * READ_WINDOW_HEIGHT);
        let plane_size = (normal_generation_count/num_threads);
        let mut modulo = normal_generation_count % num_threads;
        println!("normal slice len {}", plane_size);
        for i in 0..num_threads
        {
            let mut band_size = plane_size + 2*READ_WINDOW_HEIGHT;
            if i == 0 || i == num_threads - 1
            {
                band_size -= READ_WINDOW_HEIGHT;
            }
            if modulo > 0
            {
                band_size += 1;
                modulo -= 1;
            }
            slice_lens.push(band_size);
            println!("Slice {} len: {} ",i,band_size);
        }
        slice_lens
    }
    
    pub fn print_filled_indices(&self)
    {
        for i in 0..self.board.len()
        {
            if self.board[i][0] != Biome::Empty
            {
                println!("filled row: {i}");
            }
        }
    }

    pub fn division_with_remainder(&self, base: usize, divider: usize) ->  (usize,usize) {
            let quot = base / divider;
            let rem = base % divider;
            (quot, rem)
    }
    
    pub fn shrink_to_pow2(&self,num_sys_threads: usize) -> usize{
        
        let num_threads: usize = 2;
        let mut current_power = 0;
        while num_threads.pow(current_power + 1) <= num_sys_threads
        {
            current_power +=1;
        }
        return num_threads.pow(current_power);

    }
    
}