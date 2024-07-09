 SIMPLE unstady algorithm based CFD solver in Rust

This project is a computational fluid dynamics (CFD) solver written in Rust and post-processed with matplotlib, using the SIMPLE algorithm with a collocated grid to integrate the 2D incompressible unsteady Navier-Stokes equations. This project is based on the [lectures by Dr. Sandip Mazumder](https://youtube.com/playlist?list=PLVuuXJfoPgT4gJcBAAFPW7uMwjFKB9aqT). It can solve three cases: the lid-driven cavity flow, the pipe flow with a velocity inlet/gauge pressure outlet, and the backward facing step flow.

## Usage

To run this project, you need to have Rust, Python 3 and matplotlib installed on your system. To run the solver for a specific case, use the following command:

```bash
cargo run --release -- -c <case>
```

where `<case>` can be one of `lid_driven_cavity`, `pipe_flow`, `backward_facing_step`, `backward_facing_step`. The solver uses a uniform mesh of size $n_x \times n_y$, which can be specified by the user. To get information about what variables you can change in the solver, you can use the `--help` or `-h` flag.

## Results

https://github.com/miguelggcc/CFD-Unsteady-SIMPLE-Rust/assets/100235899/dfcf78de-3a07-4a65-9f50-2fde22fa35a4

https://github.com/miguelggcc/CFD-Unsteady-SIMPLE-Rust/assets/100235899/962ec07d-1537-42e8-b764-92460e41911f



