using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace FastFluidSolver
{
    public class FluidSolver
    {
        const int MAX_ITER = 50; //maximum number of iterations for Gauss-Seidel solver
        const double TOL = 1e-5; //maximum relative error for Gauss-Seidel solver

        public double[, ,] u { get; private set; } // x component of velocity
        public double[, ,] v { get; private set; }// y component of velocity
        public double[, ,] w { get; private set; } // z component of velocity
        public double[, ,] p { get; private set; } // pressure

        private double[, ,] u_scratch; //scratch arrays for velocities
        private double[, ,] v_scratch;
        private double[, ,] w_scratch;

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

            u_scratch = new double[Nx + 1, Ny + 1, Nz + 1];
            v_scratch = new double[Nx + 1, Ny + 1, Nz + 1];
            w_scratch = new double[Nx + 1, Ny + 1, Nz + 1];
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

            u_scratch = new double[Nx + 1, Ny + 1, Nz + 1];
            v_scratch = new double[Nx + 1, Ny + 1, Nz + 1];
            w_scratch = new double[Nx + 1, Ny + 1, Nz + 1];

            verbose = old.verbose;
        }

        /****************************************************************************
         * Diffusion step. Diffuse solve diffusion equation x_t = L(x) using second 
         * order finite difference in space and backwards Euler in time
         ****************************************************************************/
        void diffuse(ref double[, ,] x)
        {
            double[, ,] x_old = new double[x.GetLength(0), x.GetLength(1), x.GetLength(2)];
            x.CopyTo(x_old, 0);

            gs_solve(1 + 6 * nu * dt / Math.Pow(h, 2), -dt * nu / Math.Pow(h, 2), x_old, ref x, 0);
        }

        /*****************************************************************************
         * Projection step. Solves Poisson equation L(p) = div(u_old) using second order finte
         * difference and updates the velocities, u = u_old - grad(p)
         ****************************************************************************/
        void project()
        {
            double[, ,] div = new double[Nx, Ny, Nz];

            // calculate div(w) using second order finite differences
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        if (omega.obstacle_cells[i, j, k] == 0) //node not inside an obstacle
                        {
                            div[i, j, k] = (u[i + 1, j, k] - u[i, j, k]) / hx +
                                   (v[i, j + 1, k] - v[i, j, k]) / hy + (w[i, j, k + 1] - w[i, j, k]) / hz;                      
                        }
                    }
                }
            }

            gs_solve(-6 / Math.Pow(h, 2), 1 / Math.Pow(h, 2), div, ref p, 1);

            //update velocity by adding calculate grad(p), calculated using second order finite difference
            //only need to add to interior points, as the velocity at boundary points has already been fixed
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        if (omega.obstacle_cells[i, j, k] == 0) //node not inside an obstacle
                        {
                            if (omega.boundary_cells[i, j ,k] == 0) //node not on boundary, use second order finite differnce
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
        void advect(ref double[] x, double[] x0, double[] velx, double[] vely,  double[] velz)
        {
        }

        /*****************************************************************************
         * Solves the banded system given by the finite difference method applied
         * to the Poisson or diffusion equation using the iterative Gauss-Seidel method.
         * @inputs
         * double a - coefficient along diagonal entry
         * double c - coefficient of all other nonzero entries
         * double[] b - right hand side
         * double[] x - reference to array in which to store solution
         * int boundary type - 0 corresponds to Dirichlet, 1 to homogeneous Neumann 
         * 
         * TO DO: add Jacobi solver, which can be run on multiple cores
         ****************************************************************************/
        void gs_solve(double a, double[] c, double[, ,] b, ref double[, ,] x, int boundary_type)
        {
            int iter = 0;
            double res = 2 * TOL;
            while (iter < MAX_ITER && res > TOL)
            {
                res = 0;

                for (int i = 0; i < Nx + 1; i++)
                {
                    for (int j = 0; j < Ny + 1; j++)
                    {
                        for (int k = 0; k < Nz + 1; k++)
                        {
                            if (omega.obstacle_cells[i, j, k] == 0) //if node not inside obstacle
                            {
                                double x_old = x[i, j, k];

                                if (omega.boundary_cells[i, j, k] == 0) //if not on boundary, second order finite difference
                                {
                                    x[i, j, k] = (b[i, j, k] - (c[0] * x[i, j, k - 1] + c[1] * x[i, j - 1, k] + c[2] * x[i - 1, j, k] +
                                            c[3] * x[i + 1, j, k] + c[4] * x[i, j + 1, k] + c[5] * x[i, j, k + 1])) / a;
                                }
                                /*else if (boundary_type == 1)//if on boundary and homogeneous Neumann boundary conditions
                                {
                                    int nx = omega.boundary_normal_x[cell_index(i, j, k, N)];
                                    int ny = omega.boundary_normal_y[cell_index(i, j, k, N)];
                                    int nz = omega.boundary_normal_z[cell_index(i, j, k, N)];

                                    x[cell_index(i, j, k, N)] = (b[cell_index(i, j, k, N)] - c * (Math.Abs(1 + nx) * x[cell_index(i - 1, j, k, N)] +
                                            Math.Abs(1 - nx) * x[cell_index(i + 1, j, k, N)] + Math.Abs(1 + ny) * x[cell_index(i, j - 1, k, N)] +
                                            Math.Abs(1 - ny) * x[cell_index(i, j + 1, k, N)] + Math.Abs(1 + nz) * x[cell_index(i, j, k - 1, N)] +
                                            Math.Abs(1 - nz) * x[cell_index(i, j, k + 1, N)])) / a;

                                }*/

                                res += Math.Pow((x_old - x[i, j, k]), 2);
                            }
                        }
                    }
                }

                res = Math.Sqrt(res) / (Nx * Ny * Nz);
                iter++;
            }

            if (verbose)
            {
                Console.WriteLine("Gauss-Seidel solver completed with residual of {0} in {1} iterations", res, iter);
            }
        }

        /*********************************************************************************
         * Applies the boundary conditions from the domain omega to the velocities
         ********************************************************************************/
        void apply_boundary_conditions()
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    for (int k = 0; k < N; k++)
                    {
                        if (omega.boundary_nodes[cell_index(i, j, k, N)] == 1)
                        {
                            u[cell_index(i, j, k, N)] = omega.boundary_u[cell_index(i, j, k, N)];
                            v[cell_index(i, j, k, N)] = omega.boundary_v[cell_index(i, j, k, N)];
                            w[cell_index(i, j, k, N)] = omega.boundary_w[cell_index(i, j, k, N)];
                        }
                    }
                }
            }
        }

        /*****************************************************************************
         * Performs a trilinear interpolation inside a cube with sides of length h
         * @inputs:
         * double x,y,z - the x,y and z coordinates of a point inside the box, with the 
         * origin at one of the corners. All these coordinates must be between 0 and h
         * double[] values - the values at the 8 corners of the cube. These must be 
         * specified in the order:
         * 1. i+1, j+1, k+1
         * 2. i+1, j+1, k
         * 3. i+1, j, k+1
         * 4. i+1, j, k
         * 5. i, j+1, k+1
         * 6. i, j+1, k
         * 7. i, j, k+1
         * 8. i, j, k
         ****************************************************************************/
        double trilinear_interpolation(double x, double y, double z, double[] values)
        {
            double f = (values[0] * (h - x) * (h - y) * (h - z) + values[1] * (h - x) * (h - y) * z +
                values[2] * (h - x) * y * (h - z) + values[3] * (h - x) * y * z +
                values[4] * x * (h - y) * (h - z) + values[5] * x * (h - y) * z +
                values[6] * x * y * (h - z) + values[7] * x * y * z) / Math.Pow(h, 3);

                return f;
        }

        /***************************************************************************
         * Takes the x, y, z indices of a cell and returns the global coordinate
         **************************************************************************/
        public static int cell_index(int x, int y, int z, int N)
        {
            return (x >= 0 && y >= 0 && z >= 0 && x < N && y < N && z < N) ? x + y * N + z * (int)Math.Pow(N, 2) : 0;
        }

        /*****************************************************************************
         * Perform a single time step
         *****************************************************************************/
        public void time_step()
        {
            apply_boundary_conditions();

            diffuse(ref u);
            diffuse(ref v);
            diffuse(ref w);

            project();

            u.CopyTo(u_scratch, 0);
            v.CopyTo(v_scratch, 0);
            w.CopyTo(w_scratch, 0);

            advect(ref u, u_scratch, u_scratch, v_scratch, w_scratch);
            advect(ref v, v_scratch, u_scratch, v_scratch, w_scratch);
            advect(ref w, w_scratch, u_scratch, v_scratch, w_scratch);

            project();
        }
    }
}
