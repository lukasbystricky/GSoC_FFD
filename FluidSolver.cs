using System;
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

        private double dt;  //time step
        private int N;      //number of points in each coordiate direction
        private double h;   //spacing in each corrdiate direction
        private double nu;  //fluid viscosity

        void initialize();
        void add_force();

        /****************************************************************************
         * Diffusion step. Diffuse solve diffusion equation x_t = L(x) using second 
         * order finite difference in space and backwards Euler in time
         ****************************************************************************/
        void diffuse(ref double[] x)
        {
            double[] x_old = x;
            gs_solve(1 + 6 * nu * dt / Math.Pow(h, 2), dt * nu / Math.Pow(h, 2), x_old, ref x);
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
                    for (int k = 0; j < N; k++)
                    {
                        div_w[cell_index(i, j, k)] = (u[cell_index(i + 1, j, k)] - u[cell_index(i - 1, j, k)] +
                            v[cell_index(i, j + 1, k)] - v[cell_index(i, j - 1, k)] + w[cell_index(i, j, k + 1)] +
                            w[cell_index(i, j, k - 1)]) / (2 * h);
                    }
                }
            }

            gs_solve(6 / Math.Pow(h, 2), 1 / Math.Pow(h, 2), div_w, ref p);

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

        void advect();
        void export_vtk();

        /*****************************************************************************
         * Solves the banded system given by the finite difference method applied
         * to the Poisson or diffusion equation using the iterative Gauss-Seidel method.
         * @inputs
         * double a - coefficient along diagonal entry
         * double c - coefficient of all other nonzero entries
         * double[] b - right hand side
         * double[] x - reference to array in which to store solution
         * 
         * TO DO: add boundary conditions
         ****************************************************************************/
        void gs_solve(double a, double c, double[] b, ref double[] x)
        {

            int iter = 0;
            double err = 2*TOL;
            while (iter < MAX_ITER && err > TOL)
            {
                err = 0;
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < N; k++)
                        {
                            double x_old = x[cell_index(i, j, k)];
                            x[cell_index(i, j, k)] = (b[cell_index(i, j, k)] - c * (x[cell_index(i - 1, j, k)] +
                                x[cell_index(i + 1, j, k)] + x[cell_index(i, j - 1, k)] + x[cell_index(i, j + 1, k)] +
                                x[cell_index(i, j, k - 1)] + x[cell_index(i, j, k + 1)])) / a;

                            err = err + Math.Pow(x_old - x[cell_index(i, j, k)], 2)/Math.Pow(N,3);
                        }
                    }
                }

                err = Math.Sqrt(err);
            }
            
        }

        /***************************************************************************
         * Takes the x, y, z indices of a cell and returns the global coordinate
         **************************************************************************/
        int cell_index(int x, int y, int z)
        {
            return x + y * N + z * N ^ 2;
        }

    }
}
