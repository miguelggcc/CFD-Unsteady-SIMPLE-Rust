use crate::tools::{ix, tridiagonal_solver};

use super::{LidDrivenCavity, Links};

impl LidDrivenCavity {
    pub fn solver_correction(
        &self,
        x: &mut [f64],
        a_0: &[f64],
        links: &[Links],
        sources: &[f64],
        iter: usize,
        epsilon: f64,
    ) {
        let n = self.nx;
        let mut diagonal = vec![0.0; self.nx.max(self.ny)];
        let mut ax = vec![0.0; self.nx.max(self.ny)];
        let mut cx = vec![0.0; self.nx.max(self.ny) - 1];
        let mut rhs = vec![0.0; self.nx.max(self.ny)];

        for it in 0..iter {
            let mut res = 0.0;
            //Bottom left corner
            let j = 0;
            let i = 0;

            let l = &links[ix(i, j, n)];
            diagonal[i] = a_0[ix(i, j, n)];
            cx[i] = l.a_e;
            rhs[i] = -l.a_n * x[ix(i, j + 1, n)] + sources[ix(i, j, n)];
            rhs[i] += -l.a_e * x[ix(i + 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //Bottom wall
            let j = 0;

            for i in 1..self.nx - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[i] = a_0[ix(i, j, n)];
                ax[i] = l.a_w;
                cx[i] = l.a_e;
                rhs[i] = -l.a_n * x[ix(i, j + 1, n)] + sources[ix(i, j, n)];
                rhs[i] += -l.a_e * x[ix(i + 1, j, n)]
                    - l.a_w * x[ix(i - 1, j, n)]
                    - a_0[ix(i, j, n)] * x[ix(i, j, n)];
            }

            //Bottom right corner
            let j = 0;
            let i = self.nx - 1;

            let l = &links[ix(i, j, n)];
            diagonal[i] = a_0[ix(i, j, n)];
            ax[i] = l.a_w;
            rhs[i] = -l.a_n * x[ix(i, j + 1, n)] + sources[ix(i, j, n)];
            rhs[i] += -l.a_w * x[ix(i - 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            tridiagonal_solver(&ax, &diagonal, &mut cx, &mut rhs, self.nx);

            replace_row(0, x, &rhs, self.nx);

            for j in 1..self.ny - 1 {
                //Left wall
                let i = 0;

                let l = &links[ix(i, j, n)];
                diagonal[i] = a_0[ix(i, j, n)];
                cx[i] = l.a_e;
                rhs[i] =
                    -l.a_n * x[ix(i, j + 1, n)] - l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
                rhs[i] += -l.a_e * x[ix(i + 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

                //Interior nodes

                for i in 1..self.nx - 1 {
                    let l = &links[ix(i, j, n)];
                    diagonal[i] = a_0[ix(i, j, n)];
                    ax[i] = l.a_w;
                    cx[i] = l.a_e;
                    rhs[i] = -l.a_n * x[ix(i, j + 1, n)] - l.a_s * x[ix(i, j - 1, n)]
                        + sources[ix(i, j, n)];
                    rhs[i] += -l.a_e * x[ix(i + 1, j, n)]
                        - l.a_w * x[ix(i - 1, j, n)]
                        - a_0[ix(i, j, n)] * x[ix(i, j, n)];
                }

                //Right wall
                let i = self.nx - 1;

                let l = &links[ix(i, j, n)];
                diagonal[i] = a_0[ix(i, j, n)];
                ax[i] = l.a_w;
                rhs[i] =
                    -l.a_n * x[ix(i, j + 1, n)] - l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
                rhs[i] += -l.a_w * x[ix(i - 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

                tridiagonal_solver(&ax, &diagonal, &mut cx, &mut rhs, self.nx);

                replace_row(j, x, &rhs, self.nx);
            }

            //Top left corner
            let j = self.ny - 1;
            let i = 0;

            let l = &links[ix(i, j, n)];
            diagonal[i] = a_0[ix(i, j, n)];
            cx[i] = l.a_e;
            rhs[i] = -l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
            rhs[i] += -l.a_e * x[ix(i + 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //Top wall
            let j = self.ny - 1;

            for i in 1..self.nx - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[i] = a_0[ix(i, j, n)];
                ax[i] = l.a_w;
                cx[i] = l.a_e;
                rhs[i] = -l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
                rhs[i] += -l.a_e * x[ix(i + 1, j, n)]
                    - l.a_w * x[ix(i - 1, j, n)]
                    - a_0[ix(i, j, n)] * x[ix(i, j, n)];
            }

            //Top right corner
            let j = self.ny - 1;
            let i = self.nx - 1;

            let l = &links[ix(i, j, n)];
            diagonal[i] = a_0[ix(i, j, n)];
            ax[i] = l.a_w;
            rhs[i] = -l.a_s * x[ix(i, j - 1, n)] + sources[ix(i, j, n)];
            rhs[i] += -l.a_w * x[ix(i - 1, j, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            tridiagonal_solver(&ax, &diagonal, &mut cx, &mut rhs, self.nx);

            replace_row(self.ny - 1, x, &rhs, self.nx);

            //------------------------------------------------------------------------------------------------------------------------------------------------------

            //Bottom left corner
            let j = 0;
            let i = 0;

            let l = &links[ix(i, j, n)];
            diagonal[j] = a_0[ix(i, j, n)];
            cx[j] = l.a_n;
            rhs[j] = -l.a_e * x[ix(i + 1, j, n)] + sources[ix(i, j, n)];
            rhs[j] += -l.a_n * x[ix(i, j + 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //Left wall
            let i = 0;

            for j in 1..self.ny - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[j] = a_0[ix(i, j, n)];
                ax[j] = l.a_s;
                cx[j] = l.a_n;
                rhs[j] = -l.a_e * x[ix(i + 1, j, n)] + sources[ix(i, j, n)];
                rhs[j] += -l.a_n * x[ix(i, j + 1, n)]
                    - l.a_s * x[ix(i, j - 1, n)]
                    - a_0[ix(i, j, n)] * x[ix(i, j, n)];
            }

            //Top left corner
            let j = self.ny - 1;
            let i = 0;

            let l = &links[ix(i, j, n)];
            diagonal[j] = a_0[ix(i, j, n)];
            ax[j] = l.a_s;
            rhs[j] = -l.a_e * x[ix(i + 1, j, n)] + sources[ix(i, j, n)];
            rhs[j] += -l.a_s * x[ix(i, j - 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            tridiagonal_solver(&ax, &diagonal, &mut cx, &mut rhs, self.ny);

            replace_column(0, x, &rhs, self.nx, self.ny);
            add_sol_to_res(&mut res, &rhs);

            for i in 1..self.nx - 1 {
                //Bottom wall

                let j = 0;

                let l = &links[ix(i, j, n)];
                diagonal[j] = a_0[ix(i, j, n)];
                cx[j] = l.a_n;
                rhs[j] =
                    -l.a_e * x[ix(i + 1, j, n)] - l.a_w * x[ix(i - 1, j, n)] + sources[ix(i, j, n)];
                rhs[j] += -l.a_n * x[ix(i, j + 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

                //Interior nodes

                for j in 1..self.ny - 1 {
                    let l = &links[ix(i, j, n)];
                    diagonal[j] = a_0[ix(i, j, n)];
                    ax[j] = l.a_s;
                    cx[j] = l.a_n;
                    rhs[j] = -l.a_e * x[ix(i + 1, j, n)] - l.a_w * x[ix(i - 1, j, n)]
                        + sources[ix(i, j, n)];
                    rhs[j] += -l.a_n * x[ix(i, j + 1, n)]
                        - l.a_s * x[ix(i, j - 1, n)]
                        - a_0[ix(i, j, n)] * x[ix(i, j, n)];
                }

                //Top wall
                let j = self.ny - 1;
                let l = &links[ix(i, j, n)];
                diagonal[j] = a_0[ix(i, j, n)];
                ax[j] = l.a_s;
                rhs[j] =
                    -l.a_e * x[ix(i + 1, j, n)] - l.a_w * x[ix(i - 1, j, n)] + sources[ix(i, j, n)];
                rhs[j] += -l.a_s * x[ix(i, j - 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

                tridiagonal_solver(&ax, &diagonal, &mut cx, &mut rhs, self.ny);

                replace_column(i, x, &rhs, self.nx, self.ny);
                add_sol_to_res(&mut res, &rhs);
            }

            //Bottom right corner
            let j = 0;
            let i = self.nx - 1;

            let l = &links[ix(i, j, n)];
            diagonal[j] = a_0[ix(i, j, n)];
            cx[j] = l.a_n;
            rhs[j] = -l.a_w * x[ix(i - 1, j, n)] + sources[ix(i, j, n)];
            rhs[j] += -l.a_n * x[ix(i, j + 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            //Right wall
            let i = self.nx - 1;

            for j in 1..self.ny - 1 {
                let l = &links[ix(i, j, n)];
                diagonal[j] = a_0[ix(i, j, n)];
                ax[j] = l.a_s;
                cx[j] = l.a_n;
                rhs[j] = -l.a_w * x[ix(i - 1, j, n)] + sources[ix(i, j, n)];
                rhs[j] += -l.a_n * x[ix(i, j + 1, n)]
                    - l.a_s * x[ix(i, j - 1, n)]
                    - a_0[ix(i, j, n)] * x[ix(i, j, n)];
            }

            //Top right corner
            let j = self.ny - 1;
            let i = self.nx - 1;

            let l = &links[ix(i, j, n)];
            diagonal[j] = a_0[ix(i, j, n)];
            ax[j] = l.a_s;
            rhs[j] = -l.a_w * x[ix(i - 1, j, n)] + sources[ix(i, j, n)];
            rhs[j] += -l.a_s * x[ix(i, j - 1, n)] - a_0[ix(i, j, n)] * x[ix(i, j, n)];

            tridiagonal_solver(&ax, &diagonal, &mut cx, &mut rhs, self.ny);

            replace_column(self.nx - 1, x, &rhs, self.nx, self.ny);
            add_sol_to_res(&mut res, &rhs);
            if res.sqrt() < epsilon {
                dbg!(it);
                break;
            }
        }
    }
}

#[inline(always)]
fn replace_row(row: usize, x: &mut [f64], solx: &[f64], nx: usize) {
    for i in 0..nx {
        x[ix(i, row, nx)] += solx[i];
    }
}

#[inline(always)]
fn replace_column(column: usize, x: &mut [f64], soly: &[f64], nx: usize, ny: usize) {
    for j in 0..ny {
        x[ix(column, j, nx)] += soly[j];
    }
}

#[inline(always)]
fn add_sol_to_res(res: &mut f64, sol: &[f64]) {
    for s in sol {
        *res += s * s;
    }
}
