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

            double dt = 0.05;
            double nu = 1;

            double tf = 1;
            double t = 0;
        
            double[, ,] u0 = new double[Nx + 2, Ny + 2, Nz + 2];
            double[, ,] v0 = new double[Nx + 2, Ny + 2, Nz + 2];
            double[, ,] w0 = new double[Nx + 2, Ny + 2, Nz + 2];

            CavityDomain omega = new CavityDomain(Nx + 2, Ny + 2, Nz + 2);
            FluidSolver ffd = new FluidSolver(omega, dt, nu, u0, v0, w0, true);
            PostProcessor pp = new PostProcessor(ffd, omega);

            int tstep = 0;
            pp.export_vtk(String.Concat("lid_driven_cavity_", tstep, ".vtk"), Nx, Ny, Nz);
            while (t < tf)
            {
                t += dt;
                tstep++;

                Console.WriteLine("Time t = {0}", t);               
                
                ffd.time_step();
                pp.export_vtk(String.Concat("lid_driven_cavity_", tstep, ".vtk"), Nx, Ny, Nz);           
            }
        }
    }
}
