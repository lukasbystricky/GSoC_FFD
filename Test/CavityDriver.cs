using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolver
{
    /// <summary>
    /// Driver for the lid driven cavity example
    /// </summary>
    class CavityDriver
    {
        static void Main()
        {
            int Nx = 20;
            int Ny = 20;
            int Nz = 20;

            double dt = 0.01;
            double nu = 1e-2;

            double tf = 100;
            double t = 0;
        
            double[, ,] u0 = new double[Nx + 1, Ny + 2, Nz + 2];
            double[, ,] v0 = new double[Nx + 2, Ny + 1, Nz + 2];
            double[, ,] w0 = new double[Nx + 2, Ny + 2, Nz + 1];

            double[, ,] f_x = new double[Nx + 1, Ny + 2, Nz + 2];
            double[, ,] f_y = new double[Nx + 2, Ny + 1, Nz + 2];
            double[, ,] f_z = new double[Nx + 2, Ny + 2, Nz + 1];

           double  hz = 15.0 / Nz;

            FluidSolver.solver_struct solver_prams = new FluidSolver.solver_struct();

            solver_prams.tol = 1e-5;
            solver_prams.min_iter = 1;
            solver_prams.max_iter = 30;
            solver_prams.verbose = true;
            solver_prams.backtrace_order = 2;

            CavityDomain omega = new CavityDomain(Nx + 2, Ny + 2, Nz + 2);

            FluidSolver ffd = new FluidSolver(omega, dt, nu, u0, v0, w0, solver_prams);
            PostProcessor pp = new PostProcessor(ffd, omega);
            
            int tstep = 0;
            pp.export_data_vtk(String.Concat("cavity_", tstep, ".vtk"), 0, false);
            pp.export_geometry_vtk("cavity_geometry.vtk", 0);

            while (t < tf)
            {
                t += dt;
                tstep++;

                Console.WriteLine("Time t = {0}", t);   

                ffd.time_step(f_x, f_y, f_z);
                pp.export_data_vtk(String.Concat("cavity_", tstep, ".vtk"), t, false);
            }
        }
    }
}
