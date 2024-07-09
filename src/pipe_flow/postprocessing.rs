use crate::{plotter::Plot, tools::vector_delta};

use super::PipeFlow;

impl PipeFlow {
    pub fn plot(&self, plot: &mut Plot) {
        let x = vector_delta(0.0, self.width / (self.nx - 1) as f64, self.nx);
        let y = vector_delta(0.0, self.height / (self.ny - 1) as f64, self.ny);

        let vel: Vec<f64> = self
            .u
            .iter()
            .zip(&self.v)
            .map(|(u, v)| (u * u + v * v).sqrt())
            .collect();
        plot.pcolormesh(
            &x,
            &y,
            &vel,
            "plasma",
            &format!("Velocity magnitude (Re = {:.0})", self.re),
        );
        plot.xlabel("x");
        plot.ylabel("y");
    }

    pub fn update_frame(&self, plot: &mut Plot) {
        let vel: Vec<f64> = self
            .u
            .iter()
            .zip(&self.v)
            .map(|(u, v)| (u * u + v * v).sqrt())
            .collect();

        plot.update_pcolormesh(&vel, self.nx, self.ny)
    }
}
