mod plotter;

use std::time::Instant;

use clap::{arg, command, value_parser};
use unsteady_cfd_rust::{
    plotter::{Animation, Env, Plot},
    Cases,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let commands = command!()
        .args(&[
            arg!(-c --case <CASE>)
                .help("What case to solve")
                .value_parser([
                    "backward_facing_step",
                    "pipe_flow",
                    "lid_driven_cavity",
                    "square_cylinder",
                ])
                .required(true),
            arg!(-r --reynolds <RE>)
                .help("Reynolds number")
                .value_parser(value_parser!(u64))
                .default_value_if("case", "lid_driven_cavity", "100")
                .default_value_if("case", "pipe_flow", "10")
                .default_value_if("case", "backward_facing_step", "100")
                .default_value_if("case", "square_cylinder", "10"),
            arg!(-i --iterations <ITERATIONS>)
                .default_value("2500")
                .help("Max number of iterations")
                .value_parser(value_parser!(u32)),
            arg!(-x --nx <NX>)
                .help("Number of cells in the x direction")
                .value_parser(value_parser!(u32))
                .default_value_if("case", "lid_driven_cavity", "80")
                .default_value_if("case", "pipe_flow", "200")
                .default_value_if("case", "backward_facing_step", "200")
                .default_value_if("case", "square_cylinder", "250"),
            arg!(-y --ny <NY>)
                .help("Number of cells in the y direction")
                .value_parser(value_parser!(u32))
                .default_value_if("case", "lid_driven_cavity", "80")
                .default_value_if("case", "pipe_flow", "40")
                .default_value_if("case", "backward_facing_step", "40")
                .default_value_if("case", "square_cylinder", "100"),
            arg!(--relax_uv <RELAX_UV>)
                .help("Relaxation of the u,v velocites")
                .value_parser(value_parser!(f64))
                .default_value_if("case", "lid_driven_cavity", "0.8")
                .default_value_if("case", "pipe_flow", "0.8")
                .default_value_if("case", "backward_facing_step", "0.8")
                .default_value_if("case", "square_cylinder", "0.8"),
            arg!(--relax_p <RELAX_P>)
                .help("Relaxation of the pressure correction")
                .value_parser(value_parser!(f64))
                .default_value_if("case", "lid_driven_cavity", "0.2")
                .default_value_if("case", "pipe_flow", "0.2")
                .default_value_if("case", "backward_facing_step", "0.2")
                .default_value_if("case", "square_cylinder", "0.2"),
        ])
        .get_matches();

    let problem = commands
        .get_one::<String>("case")
        .expect("'Case' is required and drawing will fail if it's missing");
    let nx = *commands
        .get_one::<u32>("nx")
        .expect("'nx' is required and drawing will fail if it's missing");
    let ny = *commands
        .get_one::<u32>("ny")
        .expect("'ny' is required and drawing will fail if it's missing");
    let re = *commands
        .get_one::<u64>("reynolds")
        .expect("'Reynolds' is required and drawing will fail if it's missing");
    let relax_uv = *commands
        .get_one::<f64>("relax_uv")
        .expect("'relax_uv' is required and drawing will fail if it's missing");
    let relax_p = *commands
        .get_one::<f64>("relax_p")
        .expect("'relax_p' is required and drawing will fail if it's missing");

    let now = Instant::now();

    let dt = 2.0;
    let t_max = 600.0;

    let mut case = Cases::new(
        problem,
        nx as usize,
        ny as usize,
        re as f64,
        relax_uv,
        relax_p as f64,
        dt,
    );

    let env = Env::new();
    let plot = Plot::new(&env);

    let mut animation = Animation::new(&env, plot, 24);

    let iterations = *commands
        .get_one::<u32>("iterations")
        .expect("'iterations' is required and drawing will fail if it's missing");
    let mut iter = 0;

    let mut t = 0.0;

    let epsilon = 5e-7; //Value residuals have to be under to consider the solution converged

    case.setup_frame(&mut animation.plot);
    animation.setup_animation("out.mp4");

    while t < t_max {
        println!("Calculating t: {:.2}", t);
        while iter < iterations {
            iter += 1;
            println!("Iteration {} / {} of t: {:.2}", iter, iterations, t);

            case.iterate();
            case.update_frame(&mut animation.plot);
            animation.update_frame();

            if case.has_diverged() {
                println!("Solution diverged!");
                break;
            }
            if case.has_converged(epsilon) {
                println!("Solution has converged!");
                break;
            }
        }

        iter = 0;
        t += dt;

        case.reset();
    }
    animation.finish_animation();

    let elapsed = now.elapsed();
    println!("Elapsed: {:.2?}", elapsed);

    Ok(())
}
