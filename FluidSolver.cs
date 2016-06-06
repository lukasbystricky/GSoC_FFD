using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace FastFluidSolver
{
    /************************************************************************
     * Solves the Navier-Stokes equations using the Fast Fluid Dynamics method
     * desribed by Stam in the paper "Stable Fluids". Uses a staggered grid 
     * finite difference method to solve the spatial equations and backwards
     * Euler in time.
     * 
     * u[i, j, k] located at (i*hx, (j-0.5)*hy, (k-0.5)*hz)
     * i = 0, ..., Nx, j = 0, ... Ny + 1, k = 0, ..., Nz + 1
     * u[0, j, k] and u[Nx-1, j, k] on boundary
     * u[i, 0, k], u[i, j, 0], u[i, Ny+1, k], u[i, j, Nz+1] inside obstacle cell
     * 
     * v[i, j, k] located at ((i-0.5*hx), j*hy, (k-0.5)*hz)
     * i = 0, ..., Nx + 1, j = 0, ... Ny, k = 0, ..., Nz + 1
     * v[i, 0, k] and v[i, Ny-1, k] on boundary
     * v[0, j, k], v[i, j, 0], v[Nx+1, j, k], v[i, j, Nz+1] inside obstacle cell
     * 
     * w[i, j, k] located at((i-0.5)*hx, (j-0.5)*hy, k*hz)
     * i = 0, ..., Nx + 1, j = 0, ... Ny + 1, k = 0, ..., Nz 
     * w[i, j, 0] and w[i, j, Nz-1] on boundary
     * w[0, j, k], w[i, 0, k], w[Nx+1, j, k], w[i, Ny+1, k] inside obstacle cell
     *  
     * p[i,j,k] located at (i*hx, j*hy, k*hz)
     * i = 0, ... Nx + 1, j = 0, ... Ny + 1, k = 0, ... Nz + 1
     * p at i,j,k = 0, Nx+1 are inside obstacle cells
     ************************************************************************/
    public class FluidSolver
    {
        const int MAX_ITER = 50; //maximum number of iterations for Gauss-Seidel solver
        const double TOL = 1e-5; //maximum relative error for Gauss-Seidel solver

        public double[, ,] u; // x component of velocity
        public double[, ,] v;// y component of velocity
        public double[, ,] w; // z component of velocity
        public double[, ,] p; // pressure

        private double[, ,] u_old; //scratch arrays for velocities
        private double[, ,] v_old;
        private double[, ,] w_old;
        private double[, ,] p_old;

        private double dt;  //time step
        public int Nx { get; private set; } //number of points in each x coordinate direction
        public int Ny { get; private set; } //number of points in each y coordinate direction
        public int Nz { get; private set; } //number of points in each z coordinate direction

        private double hx;   //spacing in x coordinate direction
        private double hy;   //spacing in x coordinate direction
        private double hz;   //spacing in x coordinate direction
        private double nu;  //fluid viscosity

        private bool verbose;

        Domain omega;

        void initialize() { }
        void add_force() { }

        /****************************************************************************
         * Constructor
         ****************************************************************************/
        public FluidSolver(Domain omega, double dt, double nu, double[, ,] u0, double[, ,] v0, double[, ,] w0, bool verbose)
        {
            Nx = omega.Nx;
            Ny = omega.Ny;
            Nz = omega.Nz;

            hx = omega.hx;
            hy = omega.hy;
            hz = omega.hz;

            this.dt = dt;
            this.nu = nu;

            this.omega = omega;
            this.verbose = verbose;

            u = u0;
            v = v0;
            w = w0;

            p = new double[Nx, Ny, Nz];

            u_old = new double[u.GetLength(0), u.GetLength(1), u.GetLength(2)];
            v_old = new double[v.GetLength(0), v.GetLength(1), v.GetLength(2)];
            w_old = new double[w.GetLength(0), w.GetLength(1), w.GetLength(2)];
            p_old = new double[Nx, Ny, Nz];

            Array.Copy(u, 0, u_old, 0, u.Length);
            Array.Copy(v, 0, v_old, 0, v.Length);
            Array.Copy(w, 0, w_old, 0, w.Length);
        }

        /****************************************************************************
         * Copy constructor
         ****************************************************************************/
        public FluidSolver(FluidSolver old)
        {
            Nx = old.Nx;
            Ny = old.Ny;
            Nz = old.Nz;

            hx = old.hx;
            hy = old.hy;
            hz = old.hz;

            dt = old.dt;
            nu = old.nu;

            u = old.u;
            v = old.v;
            w = old.w;
            p = old.p;

            u_old = old.u_old;
            v_old = old.v_old;
            w_old = old.w_old;
            p_old = old.p_old;

            verbose = old.verbose;
        }

        /****************************************************************************
         * Diffusion step. Diffuse solve diffusion equation x_t = L(x) using second 
         * order finite difference in space and backwards Euler in time
         ****************************************************************************/
        void diffuse(double[, ,] x_old, ref double[, ,] x_new)
        {
            double a = 1 + 2 * nu * dt * (Math.Pow(hx, -2) + Math.Pow(hy, -2) + Math.Pow(hz, -2));
            double[] c = new double[6];

            double[, ,] b = new double[x_old.GetLength(0), x_old.GetLength(1), x_old.GetLength(2)];
            Array.Copy(x_old, 0, b, 0, x_old.Length);

            c[0] = -dt * nu * Math.Pow(hz, -2);
            c[1] = -dt * nu * Math.Pow(hy, -2);
            c[2] = -dt * nu * Math.Pow(hx, -2);
            c[3] = c[2];
            c[4] = c[1];
            c[5] = c[0];

            gs_solve(a, c, b, x_old, ref x_new);
        }

        /*****************************************************************************
         * Projection step. Solves Poisson equation L(p) = div(u_old) using second order finte
         * difference and updates the velocities, u = u_old - grad(p)
         ****************************************************************************/
        void project()
        {
            double[, ,] div = new double[Nx - 1, Ny - 1, Nz - 1];

            // calculate div(w) using second order finite differences
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        if (omega.obstacle_cells[i, j, k] == 0) //node not inside an obstacle
                        {
                            div[i, j, k] = (u[i, j, k] - u[i - 1, j, k]) / hx +
                                   (v[i, j, k] - v[i, j - 1, k]) / hy + (w[i, j, k] - w[i, j, k - 1]) / hz;                      
                        }
                    }
                }
            }

            double a = -2 * (Math.Pow(hx, -2) + Math.Pow(hy, -2) + Math.Pow(hz, -2));
            double[] c = new double[6];

            c[0] = Math.Pow(hz, -2);
            c[1] = Math.Pow(hy, -2);
            c[2] = Math.Pow(hx, -2);
            c[3] = c[2];
            c[4] = c[1];
            c[5] = c[0];

            gs_solve(a, c, div, p_old, ref p);

            //update velocity by adding calculate grad(p), calculated using second order finite difference
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        if (omega.obstacle_cells[i, j, k] == 0) //node not inside an obstacle
                        {
                            //if (omega.boundary_cells[i, j ,k] == 0)
                            {
                                u[i, j, k] -= (p[i + 1, j, k] - p[i, j, k]) / hx;
                                v[i, j, k] -= (p[i, j + 1, k] - p[i, j, k]) / hy;
                                w[i, j, k] -= (p[i, j, k + 1] - p[i, j, k]) / hz;
                            }
                        }
                    }
                }
            }
        }

        /***************************************************************************
         * Advection step. Uses a first order backtrace to update x.
         * 
         * @inputs
         * double[] x - reference to quantity to advect, can be a velocity, a
         * concentration or temperature
         * double[] x0 - initial state of x before advection
         * double[] velx - x velocity (either before advection if x is a velocity.
         * or after convection otherwise)
         * double[] vely - y velocity (either before advection if x is a velocity.
         * or after convection otherwise)
         * double[] velz - z velocity (either before advection if x is a velocity.
         * or after convection otherwise)
         **************************************************************************/
        void advect(ref double[, ,] x, double[, ,] x0, double[, ,] velx, double[, ,] vely,  double[, ,] velz, int grid_type)
        {
            int Sx = x.GetLength(0);
            int Sy = x.GetLength(1);
            int Sz = x.GetLength(2);
 
            for (int i = 1; i < Sx; i++)
            {
                for (int j = 0; j < Sy; j++)
                {
                    for (int k = 0; k < Sz; k++)
                    {
                        //find coordinate of point
                        double[] coordinate = new double[3];

                        switch (grid_type)
                        {
                            case 1: //cell centred grid
                                coordinate[0] = i * hx + 0.5;
                                coordinate[1] = j * hy + 0.5;
                                coordinate[2] = k * hz + 0.5;
                                break;

                            case 2: //x velocity
                                coordinate[0] = i * hx;
                                coordinate[1] = j * hy + 0.5;
                                coordinate[2] = k * hz + 0.5;
                                break;

                            case 3: //y velocity
                                coordinate[0] = i * hx + 0.5;
                                coordinate[1] = j * hy;
                                coordinate[2] = k * hz + 0.5;
                                break;

                            case 4: //z velocity
                                coordinate[0] = i * hx + 0.5;
                                coordinate[1] = j * hy + 0.5;
                                coordinate[2] = k * hz;
                                break;
                        }

                        //check if coordinate inside obstacle
                        int idomain = Math.Min((int)Math.Floor(coordinate[0] / hx), Nx - 1);
                        int jdomain = Math.Min((int)Math.Floor(coordinate[1] / hy), Ny - 1);
                        int kdomain = Math.Min((int)Math.Floor(coordinate[2] / hz), Nz - 1);

                        if (omega.obstacle_cells[idomain, jdomain, kdomain] == 0)
                        {
                            //interpolate velocities as needed
                            double u = Utilities.trilinear_interpolation(coordinate, velx, 2, omega);
                            double v = Utilities.trilinear_interpolation(coordinate, vely, 3, omega);
                            double w = Utilities.trilinear_interpolation(coordinate, velz, 4, omega);

                            //back step by dt
                            coordinate[0] += dt * u;
                            coordinate[1] += dt * v;
                            coordinate[2] += dt * w;

                            x[i, j, k] = Utilities.trilinear_interpolation(coordinate, x0, grid_type, omega);
                        }
                    }
                }
            }
        }

        /*********************************************************************************
         * Applies the Dirichlet boundary conditions from the domain omega to the velocities
         * and homogeneuous Neumann boundary conditions to the pressure
         ********************************************************************************/
        void apply_boundary_conditions()
        {
            //loop over all cells
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        if (omega.boundary_cells[i, j, k] == 1)
                        {
                            if (omega.boundary_normal_x[i, j, k] == -1)//boundary is on the -x side of cell
                            {
                                u[i - 1, j, k] = omega.boundary_u[i - 1, j, k];
                                v[i - 1, j, k] = 2 * omega.boundary_v[i - 1, j, k] - v[i, j, k];
                                v[i - 1, j - 1, k] = 2 * omega.boundary_v[i - 1, j, k] - v[i, j - 1, k];
                                w[i - 1, j, k] = 2 * omega.boundary_w[i - 1, j, k] - w[i, j, k];
                                w[i - 1, j, k - 1] = 2 * omega.boundary_w[i - 1, j, k] - w[i, j, k - 1];
                                p[i - 1, j, k] = p[i, j, k];
                            }

                            if (omega.boundary_normal_x[i, j, k] == 1) //boundary is to the +x side of cell
                            {
                                u[i, j, k] = omega.boundary_u[i + 1, j, k];
                                v[i + 1, j, k] = 2 * omega.boundary_v[i + 1, j, k] - v[i, j, k];
                                v[i + 1, j - 1, k] = 2 * omega.boundary_v[i + 1, j, k] - v[i, j - 1, k];
                                w[i + 1, j, k] = 2 * omega.boundary_w[i + 1, j, k] - w[i, j, k];
                                w[i + 1, j, k - 1] = 2 * omega.boundary_w[i + 1, j, k] - w[i, j, k - 1];
                                p[i + 1, j, k] = p[i, j, k];
                            }

                            if (omega.boundary_normal_y[i, j, k] == -1)//boundary is on the -y side of cell
                            {
                                u[i, j - 1, k] = 2 * omega.boundary_u[i, j - 1, k] - u[i, j, k];
                                u[i - 1, j - 1, k] = 2 * omega.boundary_u[i, j - 1, k] - u[i - 1, j, k];
                                v[i, j - 1, k] = omega.boundary_v[i, j - 1, k];
                                w[i, j - 1, k] = 2 * omega.boundary_w[i - 1, j, k] - w[i, j, k];
                                w[i, j - 1, k - 1] = 2 * omega.boundary_w[i - 1, j, k] - w[i, j, k - 1];
                                p[i, j - 1, k] = p[i, j, k];
                            }

                            if (omega.boundary_normal_y[i, j, k] == 1)//boundary is on the +y side of cell
                            {
                                u[i, j + 1, k] = 2 * omega.boundary_u[i, j + 1, k] - u[i, j, k];
                                u[i - 1, j + 1, k] = 2 * omega.boundary_u[i, j + 1, k] - u[i - 1, j, k];
                                v[i, j, k] = omega.boundary_v[i, j + 1, k];
                                w[i, j + 1, k] = 2 * omega.boundary_w[i + 1, j, k] - w[i, j, k];
                                w[i, j + 1, k - 1] = 2 * omega.boundary_w[i + 1, j, k] - w[i, j, k - 1];
                                p[i, j + 1, k] = p[i, j, k];
                            }

                            if (omega.boundary_normal_z[i, j, k] == -1)//boundary is on the -z side of cell
                            {
                                u[i, j, k - 1] = 2 * omega.boundary_u[i, j, k - 1] - u[i, j, k];
                                u[i - 1, j, k - 1] = 2 * omega.boundary_u[i, j, k - 1] - u[i - 1, j, k];
                                v[i, j, k - 1] = 2 * omega.boundary_v[i, j, k - 1] - v[i, j, k];
                                v[i, j - 1, k - 1] = 2 * omega.boundary_v[i, j, k - 1] - v[i, j - 1, k];
                                w[i, j, k - 1] = omega.boundary_w[i, j, k - 1];
                                p[i, j, k - 1] = p[i, j, k];
                            }

                            if (omega.boundary_normal_z[i, j, k] == 1)//boundary is on the +z side of cell
                            {
                                u[i, j, k + 1] = 2 * omega.boundary_u[i, j, k + 1] - u[i, j, k];
                                u[i - 1, j, k + 1] = 2 * omega.boundary_u[i, j, k + 1] - u[i - 1, j, k];
                                v[i, j, k + 1] = 2 * omega.boundary_v[i, j, k + 1] - v[i, j, k];
                                v[i, j - 1, k + 1] = 2 * omega.boundary_v[i, j, k + 1] - v[i, j - 1, k];
                                w[i, j, k] = omega.boundary_w[i, j, k + 1];
                                p[i, j, k + 1] = p[i, j, k];
                            }
                        }
                    }
                }
            }
        }


        /*****************************************************************************
         * Perform a single time step
         *****************************************************************************/
        public void time_step()
        {
            diffuse(u_old, ref u);
            diffuse(v_old, ref v);
            diffuse(w_old, ref w);

            project();

            apply_boundary_conditions();

            Array.Copy(u, 0, u_old, 0, u.Length);
            Array.Copy(v, 0, v_old, 0, v.Length);
            Array.Copy(w, 0, w_old, 0, w.Length);

           /* advect(ref u, u_old, u_old, v_old, w_old, 2);
            advect(ref v, v_old, u_old, v_old, w_old, 3);
            advect(ref w, w_old, u_old, v_old, w_old, 4);

            project();

            Array.Copy(u, 0, u_old, 0, u.Length);
            Array.Copy(v, 0, v_old, 0, v.Length);
            Array.Copy(w, 0, w_old, 0, w.Length);*/
        }

        /*****************************************************************************
         * Solves the banded system given by the finite difference method applied
         * to the Poisson or diffusion equation using the iterative Gauss-Seidel method.
         * @inputs
         * double a - coefficient along diagonal entry
         * double[] c[6] - coefficient of all other nonzero entries in order from left 
         * to right (-k, -j, -i, +i, +j, +k)
         * double[, ,] b - right hand side
         * double[, ,] x0 - initial guess
         * out double[, ,] x1 - solution
         * 
         * TO DO: add Jacobi solver, which can be run on multiple cores
         ****************************************************************************/
        void gs_solve(double a, double[] c, double[, ,] b, double[, ,] x0, ref double[, ,] x1)
        {
            int Sx = x0.GetLength(0);
            int Sy = x0.GetLength(1);
            int Sz = x0.GetLength(2);

            x1 = new double[Sx, Sy, Sz];

            int iter = 0;
            double res = 2 * TOL;

            apply_boundary_conditions();

            while (iter < MAX_ITER && res > TOL)
            {

                for (int i = 1; i < Sx - 1; i++)
                {
                    for (int j = 1; j < Sy - 1; j++)
                    {
                        for (int k = 1; k < Sz - 1; k++)
                        {
                            if (omega.obstacle_cells[i, j, k] == 0) //if cell not an obstacle
                            {
                                x1[i, j, k] = (b[i, j, k] - (c[0] * x0[i, j, k - 1] + c[1] * x0[i, j - 1, k] + c[2] * x0[i - 1, j, k] +
                                            c[3] * x0[i + 1, j, k] + c[4] * x0[i, j + 1, k] + c[5] * x0[i, j, k + 1])) / a;
                            }
                        }
                    }
                }

                apply_boundary_conditions();

                res = Utilities.compute_L2_difference(x0, x1);
                iter++;

                Array.Copy(x1, 0, x0, 0, x1.Length);
            }

            if (verbose)
            {
                Console.WriteLine("Gauss-Seidel solver completed with residual of {0} in {1} iterations", res, iter);
            }
        }
    }
}
