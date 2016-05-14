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
        const int MAX_ITER = 10; //maximum number of iterations for Jacobi solver
        const double TOL = 1e-8; //maximum relative error for Jacobi solver

        private double[] u;
        private double[] v;
        private double[] w;
        private double[] p;

        private double dt;
        private int N;
        private double h;
        private double nu;

        void initialize();
        void add_force();

        /****************************************************************************
         * Diffusion step. Diffuse solve diffusion equation for velocity, u_t = L u
         ****************************************************************************/
        void diffuse(ref double[] x)
        {
            double[] x_old = x;
            jacobi_solve(1 + 6 * nu * dt / Math.Pow(h, 2), dt * nu / Math.Pow(h, 2), x_old, ref x);
        }

        void project();
        void advect();
        void export_vtk();
        /*****************************************************************************
         * Solves the banded system given by the finite difference method applied
         * to the Poisson or diffusion equation using the iterative Jacobi method.
         * @inputs
         * double a - coefficient along diagonal entry
         * double c - coefficient of all other nonzero entries
         * double[] b - right hand side
         * double[] x - reference to array in which to store solution
         * 
         * TO DO: add boundary conditions for interior and exterior boundaries
         ****************************************************************************/
        void jacobi_solve(double a, double c, double[] b, ref double[] x)
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
