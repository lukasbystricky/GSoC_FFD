using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

/*
 * ConvergenceTest.cs
 * Copyright 2016 Lukas Bystricky <lb13f@my.fsu.edu>
 *
 * This work is licensed under the GNU GPL license version 2 or later.
 */
 
namespace FastFluidSolver
{
    /// <summary>
    /// Runs a convergence test based on the exact solution given by Ethier and Steinman.
    /// </summary>
    class Convergence
    {
        static void Main()
        {
            const double EPS = 1e-8;

            // This test varies dt and viscosity while keeping N fixed          
            int N = 64;
            double[] dt = new double[] {  1.0/50, 1.0/100, 1.0/200, 1.0/400};
            double[] nu = new double[] { 1, 0.1, 0.01, 0.001 };

            String fname = "convergence_h64_t.txt";                   

            // Constants needed for exact solution
            double a = 1.25;
            double d = 2.25;

            double tf = 1.0 / 10;

            // Create structure for solver parameters
            FluidSolver.solver_struct solver_prams = new FluidSolver.solver_struct();

            solver_prams.tol = 1e-5;
            solver_prams.min_iter = 1;
            solver_prams.max_iter = 20;
            solver_prams.verbose = false;
            solver_prams.backtrace_order = 2;

            using (StreamWriter sw = new StreamWriter(fname))
            {
                // Loop over viscosities
                for (int r = 0; r < nu.Length; r++)
                {
                    sw.WriteLine("nu = {0}", nu[r]);

                    // Loop over time step sizes
                    for (int n = 0; n < dt.Length; n++)
                    {
                        Console.WriteLine("Starting simulation for dt = {0} and nu = {1}", dt[n], nu[r]);
                        
                        sw.WriteLine("N dt  err_l2_u  err_l2_v    err_l2_w    err_l2_p    " +
                            "err_inf_u   err_inf_v   err_inf_w   err_inf_p");

                        double t = 0;
                        int Nx, Ny, Nz;
                        Nx = Ny = Nz = N;

                        double err_l2_u, err_l2_v, err_l2_w, err_l2_p;
                        double err_inf_u, err_inf_v, err_inf_w, err_inf_p;

                        double hx = 1.0 / Nx;
                        double hy = 1.0 / Ny;
                        double hz = 1.0 / Nz;

                        /*****************************************************************************************
                        * Set up initial conditions
                         *****************************************************************************************/
                        double[, ,] u0 = new double[Nx + 1, Ny + 2, Nz + 2];
                        double[, ,] v0 = new double[Nx + 2, Ny + 1, Nz + 2];
                        double[, ,] w0 = new double[Nx + 2, Ny + 2, Nz + 1];

                        double[, ,] f_x = new double[Nx + 1, Ny + 2, Nz + 2];
                        double[, ,] f_y = new double[Nx + 2, Ny + 1, Nz + 2];
                        double[, ,] f_z = new double[Nx + 2, Ny + 2, Nz + 1];

                        for (int i = 0; i < Nx + 1; i++)
                        {
                            for (int j = 0; j < Ny + 2; j++)
                            {
                                for (int k = 0; k < Nz + 2; k++)
                                {
                                    double x = i * hx;
                                    double y = (j - 0.5) * hy;
                                    double z = (k - 0.5) * hz;

                                    u0[i, j, k] = -a * (Math.Exp(a * x) * Math.Sin(a * y + d * z) +
                                        Math.Exp(a * z) * Math.Cos(a * x + d * y));
                                }
                            }
                        }

                        for (int i = 0; i < Nx + 2; i++)
                        {
                            for (int j = 0; j < Ny + 1; j++)
                            {
                                for (int k = 0; k < Nz + 2; k++)
                                {
                                    double x = (i - 0.5) * hx;
                                    double y = j * hy;
                                    double z = (k - 0.5) * hz;

                                    v0[i, j, k] = -a * (Math.Exp(a * y) * Math.Sin(a * z + d * x) +
                                        Math.Exp(a * x) * Math.Cos(a * y + d * z));
                                }
                            }
                        }

                        for (int i = 0; i < Nx + 2; i++)
                        {
                            for (int j = 0; j < Ny + 2; j++)
                            {
                                for (int k = 0; k < Nz + 1; k++)
                                {
                                    double x = (i - 0.5) * hx;
                                    double y = (j - 0.5) * hy;
                                    double z = k * hz;

                                    w0[i, j, k] = -a * (Math.Exp(a * z) * Math.Sin(a * x + d * y) +
                                        Math.Exp(a * y) * Math.Cos(a * z + d * x));
                                }
                            }
                        }

                        // Create domain
                        EthierSteinmanDomain omega = new EthierSteinmanDomain(Nx + 2, Ny + 2, Nz + 2, nu[r], a, d);

                        // Create FFD solver
                        FluidSolver ffd = new FluidSolver(omega, dt[n], nu[r], u0, v0, w0, solver_prams);
                        
                        /********************************************************************************
                         * Run time loop
                         **********************************************************************************/
                        while (Math.Abs(t - tf) > EPS)
                        {
                            t += dt[n];

                            Console.WriteLine("Time t = {0}", t);

                            // Update boundary conditions and solve single time step
                            omega.update_boundary_conditions(t);
                            ffd.time_step(f_x, f_y, f_z);
                        }

                        /*************************************************************************************
                         * Caclulate errors
                         ************************************************************************************/
                        Utilities.calculate_errors(ffd, omega, out err_l2_u, out err_inf_u, t, 2);
                        Utilities.calculate_errors(ffd, omega, out err_l2_v, out err_inf_v, t, 3);
                        Utilities.calculate_errors(ffd, omega, out err_l2_w, out err_inf_w, t, 4);
                        Utilities.calculate_errors(ffd, omega, out err_l2_p, out err_inf_p, t, 1);

                        Console.WriteLine("u : L2 error = {0}, L infinity error = {1}", err_l2_u, err_inf_u);
                        Console.WriteLine("v : L2 error = {0}, L infinity error = {1}", err_l2_v, err_inf_v);
                        Console.WriteLine("w : L2 error = {0}, L infinity error = {1}", err_l2_w, err_inf_w);
                        Console.WriteLine("p : L2 error = {0}, L infinity error = {1}", err_l2_p, err_inf_p);

                        sw.WriteLine("{0}  {1} {2} {3} {4} {5} {6} {7} {8}  {9}", N, dt[n], err_l2_u, err_l2_v, err_l2_w,
                            err_l2_p, err_inf_u, err_inf_v, err_inf_w, err_inf_p);

                        sw.Flush();
                    }

                    sw.WriteLine();
                }
            }
        }
    }
}
