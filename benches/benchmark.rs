use divan::{black_box, Bencher};
use std::sync::Arc;
use cudarc::driver::CudaDevice;
use cell_automaton::Game;
use cudarc::nvrtc::{CompileOptions, compile_ptx_with_opts};
use std::fs;

const THREADS: &[usize] = &[1, 2, 4, 8];
const SIZES: &[usize] = &[256, 512, 1024, 2048, 4096, 8192, 16384];

#[divan::bench(
    args = THREADS,
    consts = SIZES
)]
fn bench_cpu_scaling<const N:usize>(bencher: Bencher, thread_count: usize)
{
    bencher
        .with_inputs(|| Game::<N>::new())
        .bench_values(|mut game| {
            game.run(divan::black_box(thread_count));
        });
}

#[divan::bench(
    args = THREADS,
    consts = SIZES
)]
fn bench_gpu_scaling<const N:usize>(bencher: Bencher, thread_count: usize) {
    bencher
        // Używamy with_inputs, aby kompilacja i setup nie wchodziły w czas pomiaru
        .with_inputs(|| setup_cuda::<N>())
        .bench_values(|(dev, mut game)| {
            game.run_parallel_gpu(black_box(dev.clone()), 
            black_box(thread_count), 
            black_box((16,16)));
        });
}



fn setup_cuda<const N: usize >() -> (Arc<CudaDevice>, Game<N>) {
    let dev = CudaDevice::new(0).expect("NVIDIA card not found");
    let kernel_code = fs::read_to_string("kernels/automaton.cu").expect("File read error");
    
    let opts = CompileOptions {
        include_paths: vec!["C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v12.1\\include".to_string()],
        ..Default::default()
    };

    let ptx = compile_ptx_with_opts(kernel_code, opts).expect("PTX compile error");
    dev.load_ptx(ptx, "cell_automaton", &["biome_step", "setup_rand"]).expect("Error while loading kernels");

    (dev, Game::new())
}


fn main() {
    divan::main();
}