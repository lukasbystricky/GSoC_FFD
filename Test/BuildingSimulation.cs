using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/*
 * BuildingSimulation.cs
 * Copyright 2016 Lukas Bystricky <lb13f@my.fsu.edu>
 *
 * This work is licensed under the GNU GPL license version 2 or later.
 */
 
namespace FastFluidSolver
{
    /// <summary>
    /// Driver to simulate wind flow around 3 buildings.
    /// </summary>
    class BuildingSimulation
    {
        static void Main()
        {
            // Set mesh parameters, here we ask for Nx cells in the x direction
            // Ny cells in the y direction and Nz cells in the z direction (ignoring ghost cells)
            int Nx = 100;
            int Ny = 45;
            int Nz = 30;

            // Set time step and viscosity
            double dt = 0.1;
            double nu = 1e-5;

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

            // Create structure for solver parameters
            FluidSolver.solver_struct solver_prams = new FluidSolver.solver_struct();

            solver_prams.tol = 1e-3;
            solver_prams.min_iter = 1;
            solver_prams.max_iter = 30;
            solver_prams.verbose = true;
            solver_prams.backtrace_order = 2;

            // Create domain
            WindInflow omega = new WindInflow(Nx + 2, Ny + 2, Nz + 2, 100, 45, 30);

            // Add buildings
            omega.add_obstacle(15, 20, 15, 20, 0, 5);
            omega.add_obstacle(15, 18, 18, 20, 5, 12);

            omega.add_obstacle(18, 25, 25, 30, 0, 8);
            omega.add_obstacle(30, 40, 15, 25, 0, 10);

            // Create FFD solver
            FluidSolver ffd = new FluidSolver(omega, dt, nu, u0, v0, w0, solver_prams);

            // Create post processor and export initial conditions and geometry information
            PostProcessor pp = new PostProcessor(ffd, omega);
            
            int tstep = 0;
            pp.export_data_vtk(String.Concat("city_test_", tstep, ".vtk"), 0, false);
            pp.export_geometry_vtk("city_test_geometry.vtk", 0);

            // Run time loop
            while (t < tf)
            {
                t += dt;
                tstep++;

                Console.WriteLine("Time t = {0}", t);   

                // Solve single time step and export results
                ffd.time_step(f_x, f_y, f_z);
                pp.export_data_vtk(String.Concat("city_test_", tstep, ".vtk"), t, false);
            }
        }
    }
}
