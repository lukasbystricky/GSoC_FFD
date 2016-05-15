using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace FastFluidSolver
{
    class FluidSolver
    {
        const int MAX_ITER = 10; //maximum number of iterations for Gauss-Seidel solver
        const double TOL = 1e-8; //maximum relative error for Gauss-Seidel solver

        private double[] u; // x component of velocity
        private double[] v; // y component of velocity
        private double[] w; // z component of velocity
        private double[] p; // pressure

        private double[] u_sratch; //scratch arrays for velocities
        private double[] v_sratch;
        private double[] w_sratch;

        private double dt;  //time step
        private static int N;      //number of points in each coordiate direction
        private double h;   //spacing in each corrdiate direction
        private double nu;  //fluid viscosity

        Domain omega;

        void initialize() { }
        void add_force() { }

        /****************************************************************************
         * Constructor
         ****************************************************************************/
        public FluidSolver(Domain omega, double dt, double nu, double[] u0, double[] v0, double[] w0)
        {
            N = omega.N;
            h = omega.h;
            this.dt = dt;
            this.nu = nu;

            this.omega = omega;

            u = u0;
            v = v0;
            w = w0;

            p = new double[(int) Math.Pow(N, 3)];
        }
        /****************************************************************************
         * Diffusion step. Diffuse solve diffusion equation x_t = L(x) using second 
         * order finite difference in space and backwards Euler in time
         ****************************************************************************/
        void diffuse(ref double[] x)
        {
            double[] x_old = x;
            apply_boundary_conditions();
            gs_solve(1 + 6 * nu * dt / Math.Pow(h, 2), dt * nu / Math.Pow(h, 2), x_old, ref x, 0);
        }

        /*****************************************************************************
         * Projection step. Solves Poisson equation L(p) = div(u_old) using second order finte
         * difference and updates the velocities, u = u_old - grad(p)
         * 
         * TO DO: add boundary conditions
         ****************************************************************************/
        void project()
        {
            double[] div_w = new double[(int) Math.Pow(N,3)];

            // calculate div(w) using second order finite differences
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    for (int k = 0; k < N; k++)
                    {
                        div_w[cell_index(i, j, k)] = (u[cell_index(i + 1, j, k)] - u[cell_index(i - 1, j, k)] +
                            v[cell_index(i, j + 1, k)] - v[cell_index(i, j - 1, k)] + w[cell_index(i, j, k + 1)] +
                            w[cell_index(i, j, k - 1)]) / (2 * h);
                    }
                }
            }

            gs_solve(6 / Math.Pow(h, 2), 1 / Math.Pow(h, 2), div_w, ref p, 1);

            //update velocity by adding calculate grad(p), calculated using second order finite difference
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    for (int k = 0; k < N; k++)
                    {
                        u[cell_index(i, j, k)] -= (p[cell_index(i + 1, j, k)] - p[cell_index(i - 1, j, k)]) / (2 * h);
                        v[cell_index(i, j, k)] -= (p[cell_index(i, j + 1, k)] - p[cell_index(i, j - 1, k)]) / (2 * h);
                        u[cell_index(i, j, k)] -= (p[cell_index(i, j, k + 1)] - p[cell_index(i, j, k - 1)]) / (2 * h);
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
         * double[] velx - x velocity (either before convection if x is a velocity.
         * or after convection otherwise)
         * double[] vely - y velocity (either before convection if x is a velocity.
         * or after convection otherwise)
         * double[] velz - z velocity (either before convection if x is a velocity.
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
         * TO DO: add boundary conditions
         ****************************************************************************/
        void gs_solve(double a, double c, double[] b, ref double[] x, int boundary_type)
        {

            int iter = 0;
            double res = 2*TOL;
            while (iter < MAX_ITER && res > TOL)
            {
                res = 0;

                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < N; k++)
                        {
                            if (omega.obstacle[cell_index(i, j, k)] == 0) //if part of fluid domain
                            {
                                double x_old = x[cell_index(i, j, k)];

                                if (omega.boundary_nodes[cell_index(i, j, k)] == 0) //if not on boundary
                                {
                                        
                                        x[cell_index(i, j, k)] = (b[cell_index(i, j, k)] - c * (x[cell_index(i - 1, j, k)] +
                                            x[cell_index(i + 1, j, k)] + x[cell_index(i, j - 1, k)] + x[cell_index(i, j + 1, k)] +
                                            x[cell_index(i, j, k - 1)] + x[cell_index(i, j, k + 1)])) / a;
                                }
                                else if (boundary_type == 1)//if on boundary and homogeneous Neumann boundary conditions
                                {
                                        int nx = omega.boundary_normal_x[cell_index(i, j, k)];
                                        int ny = omega.boundary_normal_x[cell_index(i, j, k)];
                                        int nz = omega.boundary_normal_x[cell_index(i, j, k)];

                                        x[cell_index(i, j, k)] = (b[cell_index(i, j, k)] - c * (Math.Abs(1 + nx) * x[cell_index(i - 1, j, k)] +
                                            Math.Abs(1 - nx) * x[cell_index(i + 1, j, k)] + Math.Abs(1 + ny) * x[cell_index(i, j - 1, k)] +
                                            Math.Abs(1 - ny) * x[cell_index(i, j + 1, k)] + Math.Abs(1 + nz) * x[cell_index(i, j, k - 1)] +
                                            Math.Abs(1 - nz) * x[cell_index(i, j, k + 1)])) / a;

                                       
                                } 
                                    
                                    res = res + Math.Pow(x_old - x[cell_index(i, j, k)], 2) / Math.Pow(N, 3);
                                    break;
                                }
                            }
                        }
                    }
                }

                res = Math.Sqrt(res);
            }

        /*********************************************************************************
         * Applies the boundary conditions from omega to the velocities
         ********************************************************************************/
        void apply_boundary_conditions()
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    for (int k = 0; k < N; k++)
                    {
                        if (omega.boundary_nodes[cell_index(i, j, k)] == 1)
                        {
                            u[cell_index(i, j, k)] = omega.boundary_u[cell_index(i, j, k)];
                            v[cell_index(i, j, k)] = omega.boundary_v[cell_index(i, j, k)];
                            w[cell_index(i, j, k)] = omega.boundary_w[cell_index(i, j, k)];
                        }
                    }
                }
            }
        }

        /***************************************************************************
         * Takes the x, y, z indices of a cell and returns the global coordinate
         **************************************************************************/
        public static int cell_index(int x, int y, int z)
        {
            return (x >= 0 && y >=0 && z >= 0) ? x + y * N + z * N ^ 2 : 0;
        }

        /*****************************************************************************
         * Perform a single time step
         *****************************************************************************/
        public void time_step()
        {
            diffuse(ref u);
            diffuse(ref v);
            diffuse(ref w);

            project();
        }

        /*****************************************************************************
         * Export data to a VTK file for visualization
         ****************************************************************************/
        public void export_vtk(String fname) 
        {
            using (StreamWriter sw = new StreamWriter(fname))
            {
                sw.WriteLine("# vtk DataFile Version 3.0");
                sw.WriteLine("Fast Fluid Dynamics data\n");
                sw.WriteLine("ASCII");
                sw.WriteLine("DATASET STRUCTURED_POINTS");
                sw.WriteLine("DIMENSIONS {0} {1} {2}", N, N, N);//TO DO change to accept different domain sizes
                sw.WriteLine("ORIGIN {0} {1} {2}", 0, 0, 0);
                sw.WriteLine("SPACING {0} {1} {2}", h, h, h);

                sw.WriteLine("POINT_DATA {0}", Math.Pow(N, 3));
                sw.WriteLine("VECTORS velocity double");
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < N; k++)
                        {
                            sw.WriteLine("{0} {1} {2}", u[cell_index(i, j, k)], v[cell_index(i, j, k)], w[cell_index(i, j, k)]);
                        }
                    }
                }

            }
        }

    }
}
