using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/*
 * CavityDriver.cs
 * Copyright 2016 Lukas Bystricky <lb13f@my.fsu.edu>
 *
 * This work is licensed under the GNU GPL license version 2 or later.
 */
 
namespace FastFluidSolver
{
    /// <summary>
    /// Driver for the lid driven cavity example
    /// </summary>
    class CavityDriver
    {
        static void Main()
        {
            // Set mesh parameters, here we ask for Nx cells in the x direction
            // Ny cells in the y direction and Nz cells in the z direction (ignoring ghost cells)
            int Nx = 20;
            int Ny = 20;
            int Nz = 20;

            // Set time step and viscosity
            double dt = 0.01;
            double nu = 1e-2;

            // Set simulation start and finish times
            double tf = 100;
            double t = 0;
        
            // Set initial conditions
            double[, ,] u0 = new double[Nx + 1, Ny + 2, Nz + 2];
            double[, ,] v0 = new double[Nx + 2, Ny + 1, Nz + 2];
            double[, ,] w0 = new double[Nx + 2, Ny + 2, Nz + 1];

            // Create empty arrays for body forces
            double[, ,] f_x = new double[Nx + 1, Ny + 2, Nz + 2];
            double[, ,] f_y = new double[Nx + 2, Ny + 1, Nz + 2];
            double[, ,] f_z = new double[Nx + 2, Ny + 2, Nz + 1];

            // Create structure containing solver parameters
            FluidSolver.solver_struct solver_prams = new FluidSolver.solver_struct();

            solver_prams.tol = 1e-4;
            solver_prams.min_iter = 0;
            solver_prams.max_iter = 30;
            solver_prams.verbose = true;
            solver_prams.backtrace_order = 2;

            // Create domain
            // Pass in the number of cells in x, y and z direction (now including ghost cells)
            CavityDomain omega = new CavityDomain(Nx + 2, Ny + 2, Nz + 2);

            // Create FFD solver
            FluidSolver ffd = new FluidSolver(omega, dt, nu, u0, v0, w0, solver_prams);

            // Create post processor, export intial conditions and geometry information
            PostProcessor pp = new PostProcessor(ffd, omega);
            
            int tstep = 0;
            pp.export_data_vtk(String.Concat("cavity_", tstep, ".vtk"), 0, false);
            pp.export_geometry_vtk("cavity_geometry.vtk", 0);

            // Begin time loop
            while (t < tf)
            {
                t += dt;
                tstep++;

                Console.WriteLine("Time t = {0}", t);   

                // Solve single time step and export results
                ffd.time_step(f_x, f_y, f_z);
                pp.export_data_vtk(String.Concat("cavity_", tstep, ".vtk"), t, false);
            }
        }
    }
}
