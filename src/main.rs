use cell_automaton::{Game, game, visualization};
use cudarc::nvrtc::{CompileOptions, compile_ptx_with_opts};
use cudarc::{driver::CudaDevice, nvrtc::compile_ptx};
use criterion::{criterion_group, criterion_main, Criterion};
use std::sync::Arc;
use std::fs;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    cpu_main();
    //gpu_main();
    Ok(())
}


fn gpu_main()
{
    let (dev,mut game) = setup_cuda::<1024>();
    println!("CPU simulation start...");
    game.run_parallel_gpu(dev,4, (16,16));
    println!("End of simulation.");
    visualization::visualize(&game);
}


fn cpu_main() {
    let mut game = Game::<1024>::new(); 
    println!("CPU simulation start...");
    game.run(1);
    println!("End of simulation.");
    //visualization::visualize(&game);
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

