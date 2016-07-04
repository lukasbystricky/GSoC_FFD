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
            int Nx = 16;
            int Ny = 16;
            int Nz = 16;

            double dt = 0.1;
            double nu = 1.0/100;

            double tf = 10;
            double t = 0;
        
            double[, ,] u0 = new double[Nx + 1, Ny + 2, Nz + 2];
            double[, ,] v0 = new double[Nx + 2, Ny + 1, Nz + 2];
            double[, ,] w0 = new double[Nx + 2, Ny + 2, Nz + 1];

            FluidSolver.solver_struct solver_prams = new FluidSolver.solver_struct();

            solver_prams.tol = 1e-5;
            solver_prams.min_iter = 1;
            solver_prams.max_iter = 10;
            solver_prams.verbose = true;

            CavityDomain omega = new CavityDomain(Nx + 2, Ny + 2, Nz + 2);
            FluidSolver ffd = new FluidSolver(omega, dt, nu, u0, v0, w0, solver_prams);
            PostProcessor pp = new PostProcessor(ffd, omega);

            int tstep = 0;
            pp.export_data_vtk(String.Concat("pipe_", tstep, ".vtk"), 0, false);
            pp.export_geometry_vtk(String.Concat("pipe_geometry", tstep, ".vtk"), 0);
            while (t < tf)
            {
                t += dt;
                tstep++;

                Console.WriteLine("Time t = {0}", t);               
                
                ffd.time_step();
                pp.export_data_vtk(String.Concat("pipe_", tstep, ".vtk"), t, false);           
            }
        }
    }
}
