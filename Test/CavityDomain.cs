using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/*
 * CavityDomain.cs
 * Copyright 2016 Lukas Bystricky <lb13f@my.fsu.edu>
 *
 * This work is licensed under the GNU GPL license version 2 or later.
 */
 
namespace FastFluidSolver
{
    /// <summary>
    /// Domain for the lid driven cavity
    /// </summary>
    public class CavityDomain : Domain
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="Nx">Number of cells in x direction (including ghost cells)</param>
        /// <param name="Ny">Number of cells in y direction (including ghost cells)</param>
        /// <param name="Nz">Number of cells in z direction (including ghost cells)</param>
        public CavityDomain(int Nx, int Ny, int Nz)
        {
            this.Nx = Nx;
            this.Ny = Ny;
            this.Nz = Nz;

            length_x = 1;
            length_y = 1;
            length_z = 1;

            hx = length_x / (Nx - 2);
            hy = length_y / (Ny - 2);
            hz = length_z / (Nz - 2);            

            boundary_cells = new int[Nx, Ny, Nz];
            obstacle_cells = new int[Nx, Ny, Nz];
            boundary_normal_x = new int[Nx, Ny, Nz];
            boundary_normal_y = new int[Nx, Ny, Nz];
            boundary_normal_z = new int[Nx, Ny, Nz];
            boundary_u = new double[Nx - 1, Ny, Nz];
            boundary_v = new double[Nx, Ny - 1, Nz];
            boundary_w = new double[Nx, Ny, Nz - 1];
            outflow_boundary_x = new int[Nx, Ny, Nz];
            outflow_boundary_y = new int[Nx, Ny, Nz];
            outflow_boundary_z = new int[Nx, Ny, Nz];

            //C# default values for int or double arrays are 0, so we only need to set nonzero fields

            set_ghost_flags();
            set_boundary_flags();

            /**************************************************************************************
             * u boundary conditions
             *************************************************************************************/
            for (int i = 0; i < boundary_u.GetLength(0); i++)
            {
                for (int j = 0; j < boundary_u.GetLength(1); j++)
                {
                    boundary_u[i, j, Nz - 1] = 1;
                }
            }
        }        

        /// <summary>
        /// Copy constructor
        /// </summary>
        /// <param name="old"></param>
        public CavityDomain(CavityDomain old)
        {
            Nx = old.Nx;
            Ny = old.Ny;
            Nz = old.Nz;

            hx = old.hx;
            hy = old.hy;
            hz = old.hz;

            length_x = old.length_x;
            length_y = old.length_y;
            length_z = old.length_z;

            boundary_cells = old.boundary_cells;
            obstacle_cells = old.obstacle_cells;
            boundary_normal_x = old.boundary_normal_x;
            boundary_normal_y = old.boundary_normal_y;
            boundary_normal_z = old.boundary_normal_z;
            boundary_u = old.boundary_u;
            boundary_v = old.boundary_v;
            boundary_w = old.boundary_w;

            outflow_boundary_x = old.outflow_boundary_x;
            outflow_boundary_y = old.outflow_boundary_y;
            outflow_boundary_z = old.outflow_boundary_z;
        }
    }
}
