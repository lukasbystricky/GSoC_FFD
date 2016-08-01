using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolver
{
    class Driver
    {
        static void Main()
        {
            int Nx = 50;
            int Ny = 20;
            int Nz = 10;

            double hx = 1.0 / Nx;
            double hy = 1.0 / Ny;
            double hz = 1.0 / Nz;

            double a = 1.25;
            double d = 2.25;

            double dt = 0.01;
            double nu = 0.001;

            double tf = 1;
            double t = 0;
        
            double[, ,] u0 = new double[Nx + 1, Ny + 2, Nz + 2];
            double[, ,] v0 = new double[Nx + 2, Ny + 1, Nz + 2];
            double[, ,] w0 = new double[Nx + 2, Ny + 2, Nz + 1];

            double[, ,] f_x = new double[Nx + 1, Ny + 2, Nz + 2];
            double[, ,] f_y = new double[Nx + 2, Ny + 1, Nz + 2];
            double[, ,] f_z = new double[Nx + 2, Ny + 2, Nz + 1];

            /*for (int i = 0; i < Nx + 1; i++)
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
            }*/

            FluidSolver.solver_struct solver_prams = new FluidSolver.solver_struct();

            solver_prams.tol = 1e-5;
            solver_prams.min_iter = 1;
            solver_prams.max_iter = 20;
            solver_prams.verbose = true;
            solver_prams.backtrace_order = 1;

            TestDomain omega = new TestDomain(Nx + 2, Ny + 2, Nz + 2);
            FluidSolver ffd = new FluidSolver(omega, dt, nu, u0, v0, w0, solver_prams);
            PostProcessor pp = new PostProcessor(ffd, omega);

            int tstep = 0;
            pp.export_data_vtk(String.Concat("test_", tstep, ".vtk"), 0, false);
            pp.export_geometry_vtk(String.Concat("test_geometry", tstep, ".vtk"), 0);

            //pp.export_uninterpolated_vtk(String.Concat("pipe_u_", tstep, ".vtk"), 0);

            while (t < tf)
            {
                t += dt;
                tstep++;

                Console.WriteLine("Time t = {0}", t);   
            
                //omega.update_boundary_conditions(t);
                ffd.time_step(f_x, f_y, f_z);
                pp.export_data_vtk(String.Concat("test_", tstep, ".vtk"), t, false);
               // pp.export_uninterpolated_vtk(String.Concat("pipe_u_", tstep, ".vtk"), t);
            }
        }
    }
}
