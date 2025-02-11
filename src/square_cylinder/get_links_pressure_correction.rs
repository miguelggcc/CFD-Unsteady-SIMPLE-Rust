use crate::tools::ix;

use super::SquareCylinder;

impl SquareCylinder {
    pub fn get_links_pressure_correction(&mut self) {
        let n = self.nx;

        for j in 1..self.ny - 1 {
            for i in 1..self.nx - 1 {
                if (i < self.nx_0 || i >= self.nx_1) || (j < self.ny_0 || j >= self.ny_1) {
                    //If it is not inside the square

                    let a_0 = self.a_0[ix(i, j, n)];
                    let a_0_e = self.a_0[ix(i + 1, j, n)];
                    let a_0_w = self.a_0[ix(i - 1, j, n)];
                    let a_0_n = self.a_0[ix(i, j + 1, n)];
                    let a_0_s = self.a_0[ix(i, j - 1, n)];

                    let a_pe = -0.5 * self.dy * self.dy * (1.0 / a_0_e + 1.0 / a_0);
                    let a_pw = -0.5 * self.dy * self.dy * (1.0 / a_0_w + 1.0 / a_0);
                    let a_pn = -0.5 * self.dx * self.dx * (1.0 / a_0_n + 1.0 / a_0);
                    let a_ps = -0.5 * self.dx * self.dx * (1.0 / a_0_s + 1.0 / a_0);

                    let a_p0 = -a_pe - a_pw - a_pn - a_ps;

                    self.a_p0[ix(i, j, n)] = a_p0;

                    self.plinks[ix(i, j, n)].set_links(a_pe, a_pw, a_pn, a_ps);

                    let f = &self.faces[ix(i, j, n)];

                    let s_p = -(f.u_e - f.u_w) * self.dy - (f.v_n - f.v_s) * self.dx;

                    self.source_p[ix(i, j, n)] = s_p;
                }
            }
        }

        //---------------------Square left side---------------------------
        let i = self.nx_0 - 1;
        for j in self.ny_0..self.ny_1 {
            let a_0 = self.a_0[ix(i, j, n)];
            let a_0_w = self.a_0[ix(i - 1, j, n)];
            let a_0_n = self.a_0[ix(i, j + 1, n)];
            let a_0_s = self.a_0[ix(i, j - 1, n)];

            let a_pe = 0.0;
            let a_pw = -0.5 * self.dy * self.dy * (1.0 / a_0_w + 1.0 / a_0);
            let a_pn = -0.5 * self.dx * self.dx * (1.0 / a_0_n + 1.0 / a_0);
            let a_ps = -0.5 * self.dx * self.dx * (1.0 / a_0_s + 1.0 / a_0);
            let a_p0 = -a_pw - a_pn - a_ps;

            self.a_p0[ix(i, j, n)] = a_p0;

            self.plinks[ix(i, j, n)].set_links(a_pe, a_pw, a_pn, a_ps);

            let f = &self.faces[ix(i, j, n)];

            let s_p = -(f.u_e - f.u_w) * self.dy - (f.v_n - f.v_s) * self.dx;

            self.source_p[ix(i, j, n)] = s_p;
        }

        //---------------------Square right side---------------------------

        let i = self.nx_1;
        for j in self.ny_0..self.ny_1 {
            let a_0 = self.a_0[ix(i, j, n)];
            let a_0_e = self.a_0[ix(i + 1, j, n)];
            let a_0_n = self.a_0[ix(i, j + 1, n)];
            let a_0_s = self.a_0[ix(i, j - 1, n)];

            let a_pe = -0.5 * self.dy * self.dy * (1.0 / a_0_e + 1.0 / a_0);
            let a_pw = 0.0;
            let a_pn = -0.5 * self.dx * self.dx * (1.0 / a_0_n + 1.0 / a_0);
            let a_ps = -0.5 * self.dx * self.dx * (1.0 / a_0_s + 1.0 / a_0);

            let a_p0 = -a_pe - a_pn - a_ps;

            self.a_p0[ix(i, j, n)] = a_p0;

            self.plinks[ix(i, j, n)].set_links(a_pe, a_pw, a_pn, a_ps);

            let f = &self.faces[ix(i, j, n)];

            let s_p = -(f.u_e - f.u_w) * self.dy - (f.v_n - f.v_s) * self.dx;

            self.source_p[ix(i, j, n)] = s_p;
        }

        //---------------------Square bottom side---------------------------
        let j = self.ny_0 - 1;
        for i in self.nx_0..self.nx_1 {
            let a_0 = self.a_0[ix(i, j, n)];
            let a_0_e = self.a_0[ix(i + 1, j, n)];
            let a_0_w = self.a_0[ix(i - 1, j, n)];
            let a_0_s = self.a_0[ix(i, j - 1, n)];

            let a_pe = -0.5 * self.dy * self.dy * (1.0 / a_0_e + 1.0 / a_0);
            let a_pw = -0.5 * self.dy * self.dy * (1.0 / a_0_w + 1.0 / a_0);
            let a_pn = 0.0;
            let a_ps = -0.5 * self.dx * self.dx * (1.0 / a_0_s + 1.0 / a_0);

            let a_p0 = -a_pe - a_pw - a_ps;

            self.a_p0[ix(i, j, n)] = a_p0;

            self.plinks[ix(i, j, n)].set_links(a_pe, a_pw, a_pn, a_ps);

            let f = &self.faces[ix(i, j, n)];

            let s_p = -(f.u_e - f.u_w) * self.dy - (f.v_n - f.v_s) * self.dx;

            self.source_p[ix(i, j, n)] = s_p;
        }

        //---------------------Square top side---------------------------
        let j = self.ny_1;
        for i in self.nx_0..self.nx_1 {
            let a_0 = self.a_0[ix(i, j, n)];
            let a_0_e = self.a_0[ix(i + 1, j, n)];
            let a_0_w = self.a_0[ix(i - 1, j, n)];
            let a_0_n = self.a_0[ix(i, j + 1, n)];

            let a_pe = -0.5 * self.dy * self.dy * (1.0 / a_0_e + 1.0 / a_0);
            let a_pw = -0.5 * self.dy * self.dy * (1.0 / a_0_w + 1.0 / a_0);
            let a_pn = -0.5 * self.dx * self.dx * (1.0 / a_0_n + 1.0 / a_0);
            let a_ps = 0.0;

            let a_p0 = -a_pe - a_pw - a_pn;

            self.a_p0[ix(i, j, n)] = a_p0;

            self.plinks[ix(i, j, n)].set_links(a_pe, a_pw, a_pn, a_ps);

            let f = &self.faces[ix(i, j, n)];

            let s_p = -(f.u_e - f.u_w) * self.dy - (f.v_n - f.v_s) * self.dx;

            self.source_p[ix(i, j, n)] = s_p;
        }

        //---------------------Left wall---------------------------

        let i = 0;
        for j in 1..self.ny - 1 {
            let a_0 = self.a_0[ix(i, j, n)];
            let a_0_e = self.a_0[ix(i + 1, j, n)];
            let a_0_n = self.a_0[ix(i, j + 1, n)];
            let a_0_s = self.a_0[ix(i, j - 1, n)];

            let a_pe = -0.5 * self.dy * self.dy * (1.0 / a_0_e + 1.0 / a_0);
            let a_pw = 0.0;
            let a_pn = -0.5 * self.dx * self.dx * (1.0 / a_0_n + 1.0 / a_0);
            let a_ps = -0.5 * self.dx * self.dx * (1.0 / a_0_s + 1.0 / a_0);

            let a_p0 = -a_pe - a_pn - a_ps;

            self.a_p0[ix(i, j, n)] = a_p0;

            self.plinks[ix(i, j, n)].set_links(a_pe, a_pw, a_pn, a_ps);

            let f = &self.faces[ix(i, j, n)];

            let s_p = -(f.u_e - f.u_w) * self.dy - (f.v_n - f.v_s) * self.dx;

            self.source_p[ix(i, j, n)] = s_p;
        }

        //---------------------Right wall---------------------------
        let i = self.nx - 1;
        for j in 1..self.ny - 1 {
            let a_0 = self.a_0[ix(i, j, n)];
            let a_0_w = self.a_0[ix(i - 1, j, n)];
            let a_0_n = self.a_0[ix(i, j + 1, n)];
            let a_0_s = self.a_0[ix(i, j - 1, n)];

            let a_pe = 0.0;
            let a_pw = -0.5 * self.dy * self.dy * (1.0 / a_0_w + 1.0 / a_0)
                + 0.5 * self.dy * self.dy * (1.0 / a_0);
            let a_pn = -0.5 * self.dx * self.dx * (1.0 / a_0_n + 1.0 / a_0);
            let a_ps = -0.5 * self.dx * self.dx * (1.0 / a_0_s + 1.0 / a_0);
            let a_p0 = 0.5 * self.dy * self.dy * (1.0 / a_0_w + 1.0 / a_0)
                + 0.5 * self.dy * self.dy * (1.0 / a_0)
                - a_pn
                - a_ps;

            self.a_p0[ix(i, j, n)] = a_p0;

            self.plinks[ix(i, j, n)].set_links(a_pe, a_pw, a_pn, a_ps);

            let f = &self.faces[ix(i, j, n)];

            let s_p = -(f.u_e - f.u_w) * self.dy - (f.v_n - f.v_s) * self.dx;

            self.source_p[ix(i, j, n)] = s_p;
        }
        //---------------------Top wall---------------------------
        let j = self.ny - 1;
        for i in 1..self.nx - 1 {
            let a_0 = self.a_0[ix(i, j, n)];
            let a_0_e = self.a_0[ix(i + 1, j, n)];
            let a_0_w = self.a_0[ix(i - 1, j, n)];
            let a_0_s = self.a_0[ix(i, j - 1, n)];

            let a_pe = -0.5 * self.dy * self.dy * (1.0 / a_0_e + 1.0 / a_0);
            let a_pw = -0.5 * self.dy * self.dy * (1.0 / a_0_w + 1.0 / a_0);
            let a_pn = 0.0;
            let a_ps = -0.5 * self.dx * self.dx * (1.0 / a_0_s + 1.0 / a_0);

            let a_p0 = -a_pe - a_pw - a_ps;

            self.a_p0[ix(i, j, n)] = a_p0;

            self.plinks[ix(i, j, n)].set_links(a_pe, a_pw, a_pn, a_ps);

            let f = &self.faces[ix(i, j, n)];

            let s_p = -(f.u_e - f.u_w) * self.dy - (f.v_n - f.v_s) * self.dx;

            self.source_p[ix(i, j, n)] = s_p;
        }

        //---------------------Bottom wall---------------------------
        let j = 0;
        for i in 1..self.nx - 1 {
            let a_0 = self.a_0[ix(i, j, n)];
            let a_0_e = self.a_0[ix(i + 1, j, n)];
            let a_0_w = self.a_0[ix(i - 1, j, n)];
            let a_0_n = self.a_0[ix(i, j + 1, n)];

            let a_pe = -0.5 * self.dy * self.dy * (1.0 / a_0_e + 1.0 / a_0);
            let a_pw = -0.5 * self.dy * self.dy * (1.0 / a_0_w + 1.0 / a_0);
            let a_pn = -0.5 * self.dx * self.dx * (1.0 / a_0_n + 1.0 / a_0);
            let a_ps = 0.0;

            let a_p0 = -a_pe - a_pw - a_pn;

            self.a_p0[ix(i, j, n)] = a_p0;

            self.plinks[ix(i, j, n)].set_links(a_pe, a_pw, a_pn, a_ps);

            let f = &self.faces[ix(i, j, n)];

            let s_p = -(f.u_e - f.u_w) * self.dy - (f.v_n - f.v_s) * self.dx;

            self.source_p[ix(i, j, n)] = s_p;
        }

        //Top left corner
        let j = self.ny - 1;
        let i = 0;
        let a_0 = self.a_0[ix(i, j, n)];
        let a_0_e = self.a_0[ix(i + 1, j, n)];
        let a_0_s = self.a_0[ix(i, j - 1, n)];

        let a_pe = -0.5 * self.dy * self.dy * (1.0 / a_0_e + 1.0 / a_0);
        let a_pw = 0.0;
        let a_pn = 0.0;
        let a_ps = -0.5 * self.dx * self.dx * (1.0 / a_0_s + 1.0 / a_0);

        let a_p0 = -a_pe - a_ps;

        self.a_p0[ix(i, j, n)] = a_p0;

        self.plinks[ix(i, j, n)].set_links(a_pe, a_pw, a_pn, a_ps);

        let f = &self.faces[ix(i, j, n)];

        let s_p = -(f.u_e - f.u_w) * self.dy - (f.v_n - f.v_s) * self.dx;

        self.source_p[ix(i, j, n)] = s_p;

        //Botton left corner
        let i = 0;
        let j = 0;
        let a_0 = self.a_0[ix(i, j, n)];
        let a_0_e = self.a_0[ix(i + 1, j, n)];
        let a_0_n = self.a_0[ix(i, j + 1, n)];

        let a_pe = -0.5 * self.dy * self.dy * (1.0 / a_0_e + 1.0 / a_0);
        let a_pw = 0.0;
        let a_pn = -0.5 * self.dx * self.dx * (1.0 / a_0_n + 1.0 / a_0);
        let a_ps = 0.0;

        let a_p0 = -a_pe - a_pn;

        self.a_p0[ix(i, j, n)] = a_p0;

        self.plinks[ix(i, j, n)].set_links(a_pe, a_pw, a_pn, a_ps);

        let f = &self.faces[ix(i, j, n)];

        let s_p = -(f.u_e - f.u_w) * self.dy - (f.v_n - f.v_s) * self.dx;

        self.source_p[ix(i, j, n)] = s_p;

        //Top right corner
        let i = self.nx - 1;
        let j = self.ny - 1;
        let a_0 = self.a_0[ix(i, j, n)];
        let a_0_w = self.a_0[ix(i - 1, j, n)];
        let a_0_s = self.a_0[ix(i, j - 1, n)];

        let a_pe = 0.0;
        let a_pw = -0.5 * self.dy * self.dy * (1.0 / a_0_w + 1.0 / a_0)
            + 0.5 * self.dy * self.dy * (1.0 / a_0);
        let a_pn = 0.0;
        let a_ps = -0.5 * self.dx * self.dx * (1.0 / a_0_s + 1.0 / a_0);

        let a_p0 = 0.5 * self.dy * self.dy * (1.0 / a_0_w + 1.0 / a_0)
            + 0.5 * self.dy * self.dy * (1.0 / a_0)
            - a_ps;

        self.a_p0[ix(i, j, n)] = a_p0;

        self.plinks[ix(i, j, n)].set_links(a_pe, a_pw, a_pn, a_ps);

        let f = &self.faces[ix(i, j, n)];

        let s_p = -(f.u_e - f.u_w) * self.dy - (f.v_n - f.v_s) * self.dx;

        self.source_p[ix(i, j, n)] = s_p;

        //Bottom right corner
        let i = self.nx - 1;
        let j = 0;
        let a_0 = self.a_0[ix(i, j, n)];
        let a_0_w = self.a_0[ix(i - 1, j, n)];
        let a_0_n = self.a_0[ix(i, j + 1, n)];

        let a_pe = 0.0;
        let a_pw = -0.5 * self.dy * self.dy * (1.0 / a_0_w + 1.0 / a_0)
            + 0.5 * self.dy * self.dy * (1.0 / a_0);
        let a_pn = -0.5 * self.dx * self.dx * (1.0 / a_0_n + 1.0 / a_0);
        let a_ps = 0.0;

        let a_p0 = 0.5 * self.dy * self.dy * (1.0 / a_0_w + 1.0 / a_0)
            + 0.5 * self.dy * self.dy * (1.0 / a_0)
            - a_pn;

        self.a_p0[ix(i, j, n)] = a_p0;

        self.plinks[ix(i, j, n)].set_links(a_pe, a_pw, a_pn, a_ps);

        let f = &self.faces[ix(i, j, n)];

        let s_p = -(f.u_e - f.u_w) * self.dy - (f.v_n - f.v_s) * self.dx;

        self.source_p[ix(i, j, n)] = s_p;
    }
}
