use crate::tools::ix;

use super::SquareCylinder;

impl SquareCylinder {
    pub fn get_links_momentum(&mut self) {
        let n = self.nx;
        let inv_relax_uv = 1.0 / self.relax_uv;
        let a_0_0 = self.dx * self.dy / self.dt;

        let d_e = self.nu * self.dy / (self.dx);
        let d_w = self.nu * self.dy / (self.dx);
        let d_n = self.nu * self.dx / (self.dy);
        let d_s = self.nu * self.dx / (self.dy);

        //Interior cells
        for j in 1..self.ny - 1 {
            for i in 1..self.nx - 1 {
                if (i < self.nx_0 || i >= self.nx_1) || (j < self.ny_0 || j >= self.ny_1) {
                    //If it is not inside the square
                    let faces = &self.faces[ix(i, j, n)];

                    let a_e = -(-1.0 * faces.u_e * self.dy).max(0.0) - d_e;
                    let a_w = -(faces.u_w * self.dy).max(0.0) - d_w;
                    let a_n = -(-1.0 * faces.v_n * self.dx).max(0.0) - d_n;
                    let a_s = -(faces.v_s * self.dx).max(0.0) - d_s;

                    let a_0 = (faces.u_e * self.dy).max(0.0)
                        + d_e
                        + (-faces.u_w * self.dy).max(0.0)
                        + d_w
                        + (faces.v_n * self.dx).max(0.0)
                        + d_n
                        + (-faces.v_s * self.dx).max(0.0)
                        + d_s
                        + a_0_0;

                    self.a_0[ix(i, j, n)] = a_0 * inv_relax_uv;

                    self.links[ix(i, j, n)].set_links(a_e, a_w, a_n, a_s);

                    let s_x = 0.5 * (self.p[ix(i - 1, j, n)] - self.p[ix(i + 1, j, n)]) * self.dy;
                    let s_y = 0.5 * (self.p[ix(i, j - 1, n)] - self.p[ix(i, j + 1, n)]) * self.dx;

                    self.source_x[ix(i, j, n)] = s_x
                        + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.u[ix(i, j, n)]
                        + a_0_0 * self.u_0[ix(i, j, n)];
                    self.source_y[ix(i, j, n)] = s_y
                        + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.v[ix(i, j, n)]
                        + a_0_0 * self.v_0[ix(i, j, n)];
                }
            }
        }

        //Left side square
        let i = self.nx_0 - 1;
        for j in self.ny_0..self.ny_1 {
            let faces = &self.faces[ix(i, j, n)];

            let a_e = 0.0;
            let a_w = -(faces.u_w * self.dy).max(0.0) - d_w - d_e / 3.0;
            let a_n = -(-faces.v_n * self.dx).max(0.0) - d_n;
            let a_s = -(faces.v_s * self.dx).max(0.0) - d_s;

            let a_0 = (-faces.u_w * self.dy).max(0.0)
                + d_w
                + 3.0 * d_e
                + (faces.v_n * self.dx).max(0.0)
                + d_n
                + (-faces.v_s * self.dx).max(0.0)
                + d_s
                + a_0_0;

            self.a_0[ix(i, j, n)] = a_0 * inv_relax_uv;

            self.links[ix(i, j, n)].set_links(a_e, a_w, a_n, a_s);

            let s_x = 0.5 * (self.p[ix(i - 1, j, n)] - self.p[ix(i, j, n)]) * self.dy; //(pe+p0)/2-pe
            let s_y = 0.5 * (self.p[ix(i, j - 1, n)] - self.p[ix(i, j + 1, n)]) * self.dx;

            self.source_x[ix(i, j, n)] = s_x
                + ((1.0 - self.relax_uv) * a_0 * inv_relax_uv) * self.u[ix(i, j, n)]
                + a_0_0 * self.u_0[ix(i, j, n)];
            self.source_y[ix(i, j, n)] = s_y
                + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.v[ix(i, j, n)]
                + a_0_0 * self.v_0[ix(i, j, n)];
        }

        //Right side square
        let i = self.nx_1;
        for j in self.ny_0..self.ny_1 {
            let faces = &self.faces[ix(i, j, n)];

            let a_e = -(-faces.u_e * self.dy).max(0.0) - d_e - d_w / 3.0;
            let a_w = 0.0;
            let a_n = -(-faces.v_n * self.dx).max(0.0) - d_n;
            let a_s = -(faces.v_s * self.dx).max(0.0) - d_s;

            let a_0 = (faces.u_e * self.dy).max(0.0)
                + d_e
                + (faces.v_n * self.dx).max(0.0)
                + d_n
                + 3.0 * d_w
                + (-faces.v_s * self.dx).max(0.0)
                + d_s
                + a_0_0;

            self.a_0[ix(i, j, n)] = a_0 * inv_relax_uv;

            self.links[ix(i, j, n)].set_links(a_e, a_w, a_n, a_s);

            let s_x = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i + 1, j, n)]) * self.dy; //(p0+pw)/2-pw
            let s_y = 0.5 * (self.p[ix(i, j - 1, n)] - self.p[ix(i, j + 1, n)]) * self.dx;

            self.source_x[ix(i, j, n)] = s_x
                + ((1.0 - self.relax_uv) * a_0 * inv_relax_uv) * self.u[ix(i, j, n)]
                + a_0_0 * self.u_0[ix(i, j, n)];
            self.source_y[ix(i, j, n)] = s_y
                + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.v[ix(i, j, n)]
                + a_0_0 * self.v_0[ix(i, j, n)];
        }

        //Bottom side square
        let j = self.ny_0 - 1;
        for i in self.nx_0..self.nx_1 {
            let faces = &self.faces[ix(i, j, n)];

            let a_e = -(-faces.u_e * self.dy).max(0.0) - d_e;
            let a_w = -(faces.u_w * self.dy).max(0.0) - d_w;
            let a_n = 0.0;
            let a_s = -(faces.v_s * self.dx).max(0.0) - d_s - d_n / 3.0;

            let a_0 = (faces.u_e * self.dy).max(0.0)
                + d_e
                + (-faces.u_w * self.dy).max(0.0)
                + d_w
                + (-faces.v_s * self.dx).max(0.0)
                + d_s
                + 3.0 * d_n
                + a_0_0;

            self.a_0[ix(i, j, n)] = a_0 * inv_relax_uv;

            self.links[ix(i, j, n)].set_links(a_e, a_w, a_n, a_s);

            let s_x = 0.5 * (self.p[ix(i - 1, j, n)] - self.p[ix(i + 1, j, n)]) * self.dy;
            let s_y = 0.5 * (self.p[ix(i, j - 1, n)] - self.p[ix(i, j, n)]) * self.dx;

            self.source_x[ix(i, j, n)] = s_x
                + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.u[ix(i, j, n)]
                + a_0_0 * self.u_0[ix(i, j, n)];
            self.source_y[ix(i, j, n)] = s_y
                + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.v[ix(i, j, n)]
                + a_0_0 * self.v_0[ix(i, j, n)];
        }

        //Top side square
        let j = self.ny_1;
        for i in self.nx_0..self.nx_1 {
            let faces = &self.faces[ix(i, j, n)];

            let a_e = -(-1.0 * faces.u_e * self.dy).max(0.0) - d_e;
            let a_w = -(faces.u_w * self.dy).max(0.0) - d_w;
            let a_n = -(-1.0 * faces.v_n * self.dx).max(0.0) - d_n - d_s / 3.0;
            let a_s = 0.0;

            let a_0 = (faces.u_e * self.dy).max(0.0)
                + d_e
                + (-faces.u_w * self.dy).max(0.0)
                + d_w
                + (faces.v_n * self.dx).max(0.0)
                + d_n
                + 3.0 * d_s
                + a_0_0;

            self.a_0[ix(i, j, n)] = a_0 * inv_relax_uv;

            self.links[ix(i, j, n)].set_links(a_e, a_w, a_n, a_s);

            let s_x = 0.5 * (self.p[ix(i - 1, j, n)] - self.p[ix(i + 1, j, n)]) * self.dy;
            let s_y = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i, j + 1, n)]) * self.dx;

            self.source_x[ix(i, j, n)] = s_x
                + ((1.0 - self.relax_uv) * a_0 * inv_relax_uv) * self.u[ix(i, j, n)]
                + a_0_0 * self.u_0[ix(i, j, n)];
            self.source_y[ix(i, j, n)] = s_y
                + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.v[ix(i, j, n)]
                + a_0_0 * self.v_0[ix(i, j, n)];
        }

        //Left wall
        let i = 0;
        for j in 1..self.ny - 1 {
            let faces = &self.faces[ix(i, j, n)];

            let a_e = -(-faces.u_e * self.dy).max(0.0) - d_e - d_w / 3.0;
            let a_w = 0.0;
            let a_n = -(-faces.v_n * self.dx).max(0.0) - d_n;
            let a_s = -(faces.v_s * self.dx).max(0.0) - d_s;

            let a_0 = (faces.u_e * self.dy).max(0.0)
                + d_e
                + (faces.v_n * self.dx).max(0.0)
                + d_n
                + 3.0 * d_w
                + (-faces.v_s * self.dx).max(0.0)
                + d_s
                + a_0_0;

            self.a_0[ix(i, j, n)] = a_0 * inv_relax_uv;

            self.links[ix(i, j, n)].set_links(a_e, a_w, a_n, a_s);

            let s_x = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i + 1, j, n)]) * self.dy; //(p0+pw)/2-pw
            let s_y = 0.5 * (self.p[ix(i, j - 1, n)] - self.p[ix(i, j + 1, n)]) * self.dx;

            self.source_x[ix(i, j, n)] = s_x
                + self.u_in * self.u_in * self.dy
                + d_w * 8.0 / 3.0 * self.u_in
                + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.u[ix(i, j, n)]
                + a_0_0 * self.u_0[ix(i, j, n)];
            self.source_y[ix(i, j, n)] = s_y
                + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.v[ix(i, j, n)]
                + a_0_0 * self.v_0[ix(i, j, n)];
        }

        //Right wall
        let i = self.nx - 1;
        for j in 1..self.ny - 1 {
            let faces = &self.faces[ix(i, j, n)];

            let a_e = 0.0;
            let a_w = -(faces.u_w * self.dy).max(0.0) - d_w;
            let a_n = -(-faces.v_n * self.dx).max(0.0) - d_n;
            let a_s = -(faces.v_s * self.dx).max(0.0) - d_s;

            let a_0 = (-faces.u_w * self.dy).max(0.0)
                + d_w
                + faces.u_e * self.dy
                + (faces.v_n * self.dx).max(0.0)
                + d_n
                + (-faces.v_s * self.dx).max(0.0)
                + d_s
                + a_0_0;

            self.a_0[ix(i, j, n)] = a_0 * inv_relax_uv;

            self.links[ix(i, j, n)].set_links(a_e, a_w, a_n, a_s);

            let s_x =
                (0.5 * (self.p[ix(i - 1, j, n)] + self.p[ix(i, j, n)]) - self.p_out) * self.dy;
            let s_y = 0.5 * (self.p[ix(i, j - 1, n)] - self.p[ix(i, j + 1, n)]) * self.dx;

            self.source_x[ix(i, j, n)] = s_x
                + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.u[ix(i, j, n)]
                + a_0_0 * self.u_0[ix(i, j, n)];
            self.source_y[ix(i, j, n)] = s_y
                + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.v[ix(i, j, n)]
                + a_0_0 * self.v_0[ix(i, j, n)];
        }

        //Top wall
        let j = self.ny - 1;
        for i in 1..self.nx - 1 {
            let faces = &self.faces[ix(i, j, n)];

            let a_e = -(-faces.u_e * self.dy).max(0.0) - d_e;
            let a_w = -(faces.u_w * self.dy).max(0.0) - d_w;
            let a_n = 0.0;
            let a_s = -(faces.v_s * self.dx).max(0.0) - d_s;

            let a_0 = (faces.u_e * self.dy).max(0.0)
                + d_e
                + (-faces.u_w * self.dy).max(0.0)
                + d_w
                + (-faces.v_s * self.dx).max(0.0)
                + d_s
                + a_0_0;

            self.a_0[ix(i, j, n)] = a_0 * inv_relax_uv;

            self.links[ix(i, j, n)].set_links(a_e, a_w, a_n, a_s);

            let s_x = 0.5 * (self.p[ix(i - 1, j, n)] - self.p[ix(i + 1, j, n)]) * self.dy;
            let s_y = 0.5 * (self.p[ix(i, j - 1, n)] - self.p[ix(i, j, n)]) * self.dx;

            self.source_x[ix(i, j, n)] = s_x
                + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.u[ix(i, j, n)]
                + a_0_0 * self.u_0[ix(i, j, n)];
            self.source_y[ix(i, j, n)] = s_y
                + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.v[ix(i, j, n)]
                + a_0_0 * self.v_0[ix(i, j, n)];
        }

        //Bottom wall
        let j = 0;
        for i in 1..self.nx - 1 {
            let faces = &self.faces[ix(i, j, n)];

            let a_e = -(-1.0 * faces.u_e * self.dy).max(0.0) - d_e;
            let a_w = -(faces.u_w * self.dy).max(0.0) - d_w;
            let a_n = -(-1.0 * faces.v_n * self.dx).max(0.0) - d_n;
            let a_s = 0.0;

            let a_0 = (faces.u_e * self.dy).max(0.0)
                + d_e
                + (-faces.u_w * self.dy).max(0.0)
                + d_w
                + (faces.v_n * self.dx).max(0.0)
                + d_n
                + a_0_0;

            self.a_0[ix(i, j, n)] = a_0 * inv_relax_uv;

            self.links[ix(i, j, n)].set_links(a_e, a_w, a_n, a_s);

            let s_x = 0.5 * (self.p[ix(i - 1, j, n)] - self.p[ix(i + 1, j, n)]) * self.dy;
            let s_y = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i, j + 1, n)]) * self.dx;

            self.source_x[ix(i, j, n)] = s_x
                + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.u[ix(i, j, n)]
                + a_0_0 * self.u_0[ix(i, j, n)];
            self.source_y[ix(i, j, n)] = s_y
                + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.v[ix(i, j, n)]
                + a_0_0 * self.v_0[ix(i, j, n)];
        }

        //Top-left corner
        let j = self.ny - 1;
        let i = 0;
        let faces = &self.faces[ix(i, j, n)];

        let a_e = -(-1.0 * faces.u_e * self.dy).max(0.0) - d_e - d_w / 3.0;
        let a_w = 0.0;
        let a_n = 0.0;
        let a_s = -(faces.v_s * self.dx).max(0.0) - d_s;

        let a_0 = (faces.u_e * self.dy).max(0.0)
            + d_e
            + (-faces.v_s * self.dx).max(0.0)
            + d_s
            + 3.0 * d_w
            + a_0_0;

        self.a_0[ix(i, j, n)] = a_0 * inv_relax_uv;

        self.links[ix(i, j, n)].set_links(a_e, a_w, a_n, a_s);

        let s_x = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i + 1, j, n)]) * self.dy;
        let s_y = 0.5 * (self.p[ix(i, j - 1, n)] - self.p[ix(i, j, n)]) * self.dx;

        self.source_x[ix(i, j, n)] = s_x
            + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.u[ix(i, j, n)]
            + self.u_in * self.u_in * self.dy
            + d_w * 8.0 / 3.0 * self.u_in
            + a_0_0 * self.u_0[ix(i, j, n)];
        self.source_y[ix(i, j, n)] = s_y
            + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.v[ix(i, j, n)]
            + a_0_0 * self.v_0[ix(i, j, n)];

        //Bottom left corner
        let j = 0;
        let i = 0;
        let faces = &self.faces[ix(i, j, n)];

        let a_e = -(-1.0 * faces.u_e * self.dy).max(0.0) - d_e - d_w / 3.0;
        let a_w = 0.0;
        let a_n = -(-1.0 * faces.v_n * self.dx).max(0.0) - d_n;
        let a_s = 0.0;

        let a_0 = (faces.u_e * self.dy).max(0.0)
            + d_e
            + (faces.v_n * self.dx).max(0.0)
            + d_n
            + 3.0 * d_w
            + a_0_0;

        self.a_0[ix(i, j, n)] = a_0 * inv_relax_uv;

        self.links[ix(i, j, n)].set_links(a_e, a_w, a_n, a_s);

        let s_x = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i + 1, j, n)]) * self.dy;
        let s_y = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i, j + 1, n)]) * self.dx;

        self.source_x[ix(i, j, n)] = s_x
            + self.u_in * self.u_in * self.dy
            + d_w * 8.0 / 3.0 * self.u_in
            + ((1.0 - self.relax_uv) * a_0 * inv_relax_uv) * self.u[ix(i, j, n)]
            + a_0_0 * self.u_0[ix(i, j, n)];
        self.source_y[ix(i, j, n)] = s_y
            + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.v[ix(i, j, n)]
            + a_0_0 * self.v_0[ix(i, j, n)];

        //Top right corner
        let j = self.ny - 1;
        let i = self.nx - 1;
        let faces = &self.faces[ix(i, j, n)];

        let a_e = 0.0;
        let a_w = -(faces.u_w * self.dy).max(0.0) - d_w;
        let a_n = 0.0;
        let a_s = -(faces.v_s * self.dx).max(0.0) - d_s;

        let a_0 = faces.u_e * self.dy
            + (-faces.u_w * self.dy).max(0.0)
            + d_w
            + (-faces.v_s * self.dx).max(0.0)
            + d_s
            + a_0_0;

        self.a_0[ix(i, j, n)] = a_0 * inv_relax_uv;

        self.links[ix(i, j, n)].set_links(a_e, a_w, a_n, a_s);

        let s_x = (0.5 * (self.p[ix(i - 1, j, n)] + self.p[ix(i, j, n)]) - self.p_out) * self.dy;
        let s_y = 0.5 * (self.p[ix(i, j - 1, n)] - self.p[ix(i, j, n)]) * self.dx;

        self.source_x[ix(i, j, n)] = s_x
            + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.u[ix(i, j, n)]
            + a_0_0 * self.u_0[ix(i, j, n)];
        self.source_y[ix(i, j, n)] = s_y
            + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.v[ix(i, j, n)]
            + a_0_0 * self.v_0[ix(i, j, n)];

        //Bottom right corner
        let j = 0;
        let i = self.nx - 1;
        let faces = &self.faces[ix(i, j, n)];

        let a_e = 0.0;
        let a_w = -(faces.u_w * self.dy).max(0.0) - d_w;
        let a_n = -(-1.0 * faces.v_n * self.dx).max(0.0) - d_n;
        let a_s = 0.0;

        let a_0 = faces.u_e * self.dy
            + (-faces.u_w * self.dy).max(0.0)
            + d_w
            + (faces.v_n * self.dx).max(0.0)
            + d_n
            + a_0_0;

        self.a_0[ix(i, j, n)] = a_0 * inv_relax_uv;

        self.links[ix(i, j, n)].set_links(a_e, a_w, a_n, a_s);

        let s_x = (0.5 * (self.p[ix(i - 1, j, n)] + self.p[ix(i, j, n)]) - self.p_out) * self.dy;
        let s_y = 0.5 * (self.p[ix(i, j, n)] - self.p[ix(i, j + 1, n)]) * self.dx;

        self.source_x[ix(i, j, n)] = s_x
            + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.u[ix(i, j, n)]
            + a_0_0 * self.u_0[ix(i, j, n)];
        self.source_y[ix(i, j, n)] = s_y
            + (1.0 - self.relax_uv) * a_0 * inv_relax_uv * self.v[ix(i, j, n)]
            + a_0_0 * self.v_0[ix(i, j, n)];
    }
}
