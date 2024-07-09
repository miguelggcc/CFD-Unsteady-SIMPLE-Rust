mod backward_facing_step;
pub mod lid_driven_cavity;
mod pipe_flow;
pub mod plotter;
mod square_cylinder;
mod tools;

use backward_facing_step::BackwardFacingStep;
use lid_driven_cavity::LidDrivenCavity;
use pipe_flow::PipeFlow;
use square_cylinder::SquareCylinder;

use crate::plotter::Plot;

pub trait Case {
    fn iterate(&mut self);
    fn reset(&mut self);
    fn has_converged(&self, epsilon: f64) -> bool;
    fn has_diverged(&self) -> bool;
    fn setup_frame(&self, plot: &mut Plot);
    fn update_frame(&self, plot: &mut Plot);
}

pub enum Cases {
    LidDrivenCavity(LidDrivenCavity),
    PipeFlow(PipeFlow),
    BackwardFacingStep(BackwardFacingStep),
    SquareCylinder(SquareCylinder),
}

impl Cases {
    pub fn new(
        case: &str,
        nx: usize,
        ny: usize,
        re: f64,
        relax_uv: f64,
        relax_p: f64,
        dt: f64,
    ) -> Self {
        match case {
            "lid_driven_cavity" => {
                Cases::LidDrivenCavity(LidDrivenCavity::new(nx, ny, re, relax_uv, relax_p, dt))
            }
            "pipe_flow" => Cases::PipeFlow(PipeFlow::new(nx, ny, re, relax_uv, relax_p, dt)),
            "backward_facing_step" => Cases::BackwardFacingStep(BackwardFacingStep::new(
                nx, ny, re, relax_uv, relax_p, dt,
            )),
            "square_cylinder" => {
                Cases::SquareCylinder(SquareCylinder::new(nx, ny, re, relax_uv, relax_p, dt))
            }
            _ => unreachable!(),
        }
    }

    pub fn iterate(&mut self) {
        match self {
            Cases::LidDrivenCavity(lid_driven_cavity) => lid_driven_cavity.iterate(),
            Cases::PipeFlow(pipe_flow) => pipe_flow.iterate(),
            Cases::BackwardFacingStep(backward_facing_step) => backward_facing_step.iterate(),
            Cases::SquareCylinder(square_cylinder) => square_cylinder.iterate(),
        }
    }

    pub fn reset(&mut self) {
        match self {
            Cases::LidDrivenCavity(lid_driven_cavity) => lid_driven_cavity.reset(),
            Cases::PipeFlow(pipe_flow) => pipe_flow.reset(),
            Cases::BackwardFacingStep(backward_facing_step) => backward_facing_step.reset(),
            Cases::SquareCylinder(square_cylinder) => square_cylinder.reset(),
        }
    }

    pub fn has_converged(&self, epsilon: f64) -> bool {
        match self {
            Cases::LidDrivenCavity(lid_driven_cavity) => lid_driven_cavity.has_converged(epsilon),
            Cases::PipeFlow(pipe_flow) => pipe_flow.has_converged(epsilon),
            Cases::BackwardFacingStep(backward_facing_step) => {
                backward_facing_step.has_converged(epsilon)
            }
            Cases::SquareCylinder(square_cylinder) => square_cylinder.has_converged(epsilon),
        }
    }

    pub fn has_diverged(&self) -> bool {
        match self {
            Cases::LidDrivenCavity(lid_driven_cavity) => lid_driven_cavity.has_diverged(),
            Cases::PipeFlow(pipe_flow) => pipe_flow.has_diverged(),
            Cases::BackwardFacingStep(backward_facing_step) => backward_facing_step.has_diverged(),
            Cases::SquareCylinder(square_cylinder) => square_cylinder.has_diverged(),
        }
    }

    pub fn setup_frame(&self, plot: &mut Plot) {
        match self {
            Cases::LidDrivenCavity(lid_driven_cavity) => lid_driven_cavity.setup_frame(plot),
            Cases::PipeFlow(pipe_flow) => pipe_flow.setup_frame(plot),
            Cases::BackwardFacingStep(backward_facing_step) => {
                backward_facing_step.setup_frame(plot)
            }
            Cases::SquareCylinder(square_cylinder) => square_cylinder.setup_frame(plot),
        }
    }

    pub fn update_frame(&self, plot: &mut Plot) {
        match self {
            Cases::LidDrivenCavity(lid_driven_cavity) => lid_driven_cavity.update_frame(plot),
            Cases::PipeFlow(pipe_flow) => pipe_flow.update_frame(plot),
            Cases::BackwardFacingStep(backward_facing_step) => {
                backward_facing_step.update_frame(plot)
            }
            Cases::SquareCylinder(square_cylinder) => square_cylinder.update_frame(plot),
        }
    }
}
