use rand::Rng;
use rand::rand_core::block;
use crate::{ Biome, COAST_BIOME_CHANCE, HIGHLANDS_BIOME_CHANCE, LOWLANDS_BIOME_CHANCE, MAP_SIZE, MOUNTAINS_BIOME_CHANCE, OCEAN_BIOME_CHANCE, READ_WINDOW_HEIGHT, READ_WINDOW_WIDTH, SAFETY_OFFSET_SIZE, SEA_BIOME_CHANCE};
use cudarc::driver::{CudaDevice, CudaFunction, CudaSlice, CudaStream, DeviceRepr, LaunchAsync, LaunchConfig};
use core::slice;
use std::collections::VecDeque;
use std::os::windows::thread;
use std::str::SplitWhitespace;
use std::sync::Arc;
//N is the mapsize
pub struct Game<const N: usize>  {
    pub board: Vec<[Biome; N]>,
    pub map_size: usize
}


#[repr(C)]
#[derive(Clone, Copy, Default, Debug)]
pub struct GpuCoord {
    pub y: i32,
    pub x: i32,
}

// Implementujemy DeviceRepr dla naszej struktury
unsafe impl DeviceRepr for GpuCoord {}


#[derive(Clone, Copy)]
pub struct LimitOffsets
{
    pub read_up_offset: usize,
    pub read_down_offset: usize
}

impl<const N: usize> Game<N> {
    
    /*Constructor initializing the board */
    pub fn new() -> Self {
        let empty_row = [Biome::Empty; N];
        let board = vec![empty_row; N];
        return Game {board: board, map_size: N};
    }    

    pub fn set_thread_count(&self) -> usize
    {
        let num_sys_threads = num_cpus::get();
        let num_threads = self.shrink_to_pow2(num_sys_threads);
        num_threads
    }


    pub fn run(&mut self, thread_count: usize){
        
        let offsets = self.calculate_offsets(thread_count);
        if thread_count == 1
        {
            Game::process_slice_worker(&mut self.board, LimitOffsets{read_down_offset:0, read_up_offset:0});
            return;
        }
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
    pub fn run_parallel_gpu(&mut self, dev: Arc<CudaDevice>, spawn_count: usize, (block_size_x,block_size_y): (u32,u32)) {
        let step_fn = dev.get_func("cell_automaton", "biome_step").expect("biome_step not found");
        let setup_fn = dev.get_func("cell_automaton", "setup_rand").expect("setup_rand not found");
        
        //calculate spawn lens
        let offsets = self.calculate_offsets(spawn_count);
        let slice_lens = self.calculate_normal_slice_lens(spawn_count);
        println!("{slice_lens:?}");
        let mut y_index = 0;
        for i in 0..spawn_count
        {
            self.board[y_index + slice_lens[i]/2][self.map_size/2] = Biome::Ocean;
            y_index += slice_lens[i]
        }

        let width = N as i32;
        let height = self.board.len() as i32;
        let num_pixels = N * self.board.len();


        let flat_data: Vec<u8> = self.board.iter().flatten().map(|&b| b as u8).collect();
        //5 bytes + 1 for states times number of pixels

        let mut rng_states = dev.alloc_zeros::<u8>(num_pixels * 48).expect("Allocation error for RNG states");

        let cfg_rng = LaunchConfig::for_num_elems(num_pixels as u32);
        unsafe { 
            setup_fn.launch(cfg_rng, (&mut rng_states, 67u32, width, height)) 
        }.expect("Could not launch setup_rand kernel");
        
        let mut d_in = dev.htod_copy(flat_data.clone()).expect("Allocation error for slice data input");
        let mut d_out = dev.htod_copy(flat_data).expect("Allocation error for slice data output");
        
        //lil trick for guarenteeing that there are always threads when map is f.e. width 20 height 20
        let grid_x = (width as u32 + block_size_x - 1) / block_size_x;
        let grid_y = (height as u32 + block_size_y - 1) / block_size_y;
        let cfg = LaunchConfig {
            grid_dim: (
                grid_x,
                grid_y,
                1,
            ),
            block_dim: (block_size_x, block_size_y, 1),
            shared_mem_bytes: 0,
        };

        let mut borders = Vec::<i32>::new();
        let mut border_y_index = 0;

        for i in 0..spawn_count-1
        {
            border_y_index += (self.map_size/spawn_count) as i32;
            //start inclusive in cuda, end exclusicve
            let zone_start = border_y_index - SAFETY_OFFSET_SIZE as i32;
            let zone_end = border_y_index + SAFETY_OFFSET_SIZE as i32;
            borders.push(zone_start); 
            borders.push(zone_end);
        }
        println!("Doctors with borders {:?}", borders);
        let mut d_borders = dev.htod_copy(borders.clone()).expect("Allocation error for slice data output");
        
        let normal_loop_size = (N/2 + 1) + *slice_lens.iter().max().unwrap() + 1;
        let safety_loop_size = (N/2 + 1) + SAFETY_OFFSET_SIZE + 1;
        for i in 0..normal_loop_size{
            unsafe {
                step_fn.clone().launch(cfg, (
                    &d_in,     
                    &mut d_out,
                    &mut rng_states,
                    &d_borders,
                    borders.len(),
                    width,
                    height,
                    READ_WINDOW_HEIGHT as i32,
                    READ_WINDOW_WIDTH as i32,
                    true, //normal phase
                )).expect("Could not launch step kernel");
            }
            std::mem::swap(&mut d_in, &mut d_out);
        }
        for i in 0..safety_loop_size{
            unsafe {
                step_fn.clone().launch(cfg, (
                    &d_in,     
                    &mut d_out,
                    &mut rng_states,
                    &d_borders,
                    borders.len(),
                    width,
                    height,
                    READ_WINDOW_HEIGHT as i32,
                    READ_WINDOW_WIDTH as i32,
                    false, //safety phase
                )).expect("Could not launch step kernel");
            }
            std::mem::swap(&mut d_in, &mut d_out);
        }
        let final_flat = dev.dtoh_sync_copy(&d_in).expect("Błąd pobierania danych z GPU");

        //map to slice
        for (row_idx, row_data) in final_flat.chunks(N).enumerate() {
            for (col_idx, &val) in row_data.iter().enumerate() {
                self.board[row_idx][col_idx] = unsafe { std::mem::transmute(val) };
            }
        }

    }


    fn count_filled_pixels(slice: &[[Biome; N]], width: usize, height: usize) -> usize {
        let mut count = 0;
        for i in 0..height {
            for j in 0..width {
                if slice[i][j] != Biome::Empty {
                    count += 1;
                }
            }
        }
        println!("Liczba wypełnionych pikseli w slice: {}", count);
        count
    }
    
    
    pub fn process_slice_worker(slice: &mut [[Biome; N]], offsets: LimitOffsets)
    {
        let mut transform_indices = VecDeque::<(usize,usize)>::from([(slice.len()/2, N/2)]);
        while !transform_indices.is_empty()
        {
            let mut next_layer = Vec::<(usize,usize)>::new();
            let indices_len = transform_indices.len();
            for _ in 0..indices_len
            {
                let current = transform_indices.pop_back().unwrap();
                let adjacent_cells = Game::get_adjacent_cell_indices(slice, current, offsets);
                for (i,j) in &adjacent_cells {
                    slice[*i][*j] = Biome::Checked;
                    next_layer.push((*i,*j));
                }
                for (i,j) in adjacent_cells
                {
                    slice[i][j] = Game::calculate_biome((i,j), slice);
                }
            }
            for n in next_layer {
                transform_indices.push_back(n);
            }
        }
    }

    pub unsafe fn setup_slices(&mut self, thread_count: usize) -> (Vec<&mut [[Biome;N]]>, Vec<&mut [[Biome;N]]>)
    {
        let normal_band_sizes = self.calculate_normal_slice_lens(thread_count);
        let mut normal_bands = Vec::<&mut [[Biome;N]]>::new();
        let mut safety_bands = Vec::<&mut [[Biome;N]]>::new();
        let mut index: usize = 0;
        for i in 0..thread_count
        {
            let mut start = index;
            let mut end = index + normal_band_sizes[i];
            //println!("Normal start: {}, end {}", start, end);
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
                //println!("Safety start: {}, end {}", start, end);
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
    // Returns indices of all neighbours within the read window
    pub fn get_neighbours(slice: &mut [[Biome; N]],index: (usize,usize)) -> Vec::<(usize,usize)>
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
                    && nj < N
                    {
                        neighbours.push((ni,nj));
                    }   
                }  
            }
        }
        neighbours
    }
    
    pub fn calculate_biome(index: (usize,usize), slice: &mut [[Biome; N]]) -> Biome
    {
        let neighbours = Game::get_neighbours(slice, index);
        let mut chances: [f64; 6] = [0.0; 6];
        let mut mult: f64;

        for (i, j) in neighbours {
            
            if slice[i][j] == Biome::Empty {
                continue;
            }
            
            mult = if (i as i32).abs() + (j as i32).abs() > 1 {0.25} else {1.0};
            match slice[i][j] {
                Biome::Ocean => {
                    for i in 0..6 {
                        chances[i] += OCEAN_BIOME_CHANCE[i] * mult;
                    }
                },
                Biome::Sea => {
                    for i in 0..6 {
                        chances[i] += SEA_BIOME_CHANCE[i] * mult;
                    }
                },
                Biome::Coast => {
                    for i in 0..6 {
                        chances[i] += COAST_BIOME_CHANCE[i] * mult;
                    }
                },
                Biome::Lowlands => {
                    for i in 0..6 {
                        chances[i] += LOWLANDS_BIOME_CHANCE[i] * mult;
                    }
                },
                Biome::Highlands => {
                    for i in 0..6 {
                        chances[i] += HIGHLANDS_BIOME_CHANCE[i] * mult;
                    }
                },
                Biome::Mountains => {
                    for i in 0..6 {
                        chances[i] += MOUNTAINS_BIOME_CHANCE[i] * mult;
                    }
                },
                _ => {}
            }
        } 

        for val in &mut chances {
            if *val < 0.0 {
                *val = 0.0;
            }
        }
        
        let mut rng = rand::rng();
        let rand_val: f64 = rng.random();
        
        let sum: f64 = chances.iter().sum();

        let mut selected = Biome::Empty;
        if sum <= 0.0 {
            let mut running_sum: f64 = 0.0; 
            for i in 0..6 {
                running_sum += 1.0 / 6.0;
                
                if rand_val <= running_sum {
                    selected = Biome::from_u8((i + 1) as u8);
                    break;
                }
            }
        }
        else {
            let mut running_sum: f64 = 0.0; 
            for i in 0..6 {
                running_sum += chances[i] / sum;
                if rand_val <= running_sum {
                    selected = Biome::from_u8((i + 1) as u8);
                    break;
                }
            }
        }
        selected
    }

    // Returns top bottom/left/right adjacent cells indices that were empty and marks them as checked
    pub fn get_adjacent_cell_indices(slice: &mut [[Biome; N]], index: (usize, usize), offset: LimitOffsets) -> Vec<(usize, usize)>
    {
        let mut adjacent_cell_indices = Vec::<(usize,usize)>::new();
        let offsets = Vec::from([(0,1), (1,0), (0,-1), (-1,0)]);
        for (oi, oj) in offsets
        {
            if let (Some(ni), Some(nj)) = 
            (index.0.checked_add_signed(oi), index.1.checked_add_signed(oj)) {
                if ni >= READ_WINDOW_HEIGHT - offset.read_up_offset 
                && ni < slice.len() - (READ_WINDOW_HEIGHT - offset.read_down_offset) 
                && nj < N {
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
        let normal_generation_count = N - (safety_thread_count * SAFETY_OFFSET_SIZE);
        let plane_size = (normal_generation_count/num_threads);
        let mut modulo = normal_generation_count % num_threads;
        //println!("normal slice len {}", plane_size);
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
            //println!("Slice {} len: {} ",i,band_size);
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
    pub fn print_board(&self) {
        for i in 0..self.board.len() {
            for j in 0..self.board[i].len() {
                let symbol = match self.board[i][j] {
                    Biome::Ocean => "O", // Ocean
                    Biome::Sea => "S", // Morze (Sea)
                    Biome::Coast => "C", // Wybrzeże (Coast)
                    Biome::Lowlands => "L", // Niziny (Lowlands)
                    Biome::Highlands => "H", // Wyżyny (Highlands)
                    Biome::Mountains => "M", // Góry (Mountains)
                    _ => ".", // Puste / Brak biomu
                };
                print!("{} ", symbol);
            }
            println!();
        }
}
    
}