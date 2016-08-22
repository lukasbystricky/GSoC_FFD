using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/*
 * WindInflow.cs
 * Copyright 2016 Lukas Bystricky <lb13f@my.fsu.edu>
 *
 * This work is licensed under the GNU GPL license version 2 or later.
 */
 
namespace FastFluidSolver
{
    /// <summary>
    /// Domain with an exponential wind profile on the inflow (x = 0), 
    /// 0 velocity on the ground (z = 0) and all other boundaries marked as outflow.
    /// </summary>
    public class WindInflow : Domain
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="Nx">Number of cells (including ghost cells) in x direction</param>
        /// <param name="Ny">Number of cells (including ghost cells) in x direction</param>
        /// <param name="Nz">Number of cells (including ghost cells) in x direction</param>
        /// <param name="length_x">Length of domain in x direction (not including ghost cells)</param>
        /// <param name="length_y">Length of domain in y direction (not including ghost cells)</param>
        /// <param name="length_z">Length of domain in z direction (not including ghost cells)</param>
        public WindInflow(int Nx, int Ny, int Nz, double length_x, 
            double length_y, double length_z)
        {
            this.Nx = Nx;
            this.Ny = Ny;
            this.Nz = Nz;

            this.length_x = length_x;
            this.length_y = length_y;
            this.length_z = length_z;

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

            add_obstacle(15, 20, 15, 20, 0, 5);
            add_obstacle(15, 18, 18, 20, 5, 12);

            add_obstacle(18, 25, 25, 30, 0, 8);
            add_obstacle(30, 40, 15, 25, 0, 10);

            //set outflow boundaries
            //x = 0 will be inflow, z = 0 will be solid ground, all others will be outflow            
            
            //x outflow
            for (int j = 1; j < Ny - 1; j++)
            {
                for (int k = 1; k < Nz - 1; k++)
                {
                    outflow_boundary_x[Nx - 2, j, k] = 1;
                }
            }

            //y outflows
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int k = 1; k < Nz - 1; k++)
                {
                    outflow_boundary_y[i, 1, k] = 1;
                    outflow_boundary_y[i, Ny - 2, k] = 1;
                }
            }

            //z outflow
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    outflow_boundary_z[i, j, Nz - 2] = 1;
                }
            }

            //x = 0, inflow
            for (int j = 0; j < boundary_u.GetLength(1); j++)
            {
                for (int k = 0; k < boundary_u.GetLength(2); k++)
                {
                    double z = (k - 0.5) * hz;
                    boundary_u[0, j, k] = 6.751 * Math.Pow(Math.Max((z / 7.5), 0), 1.0 / 5.25);
                }
            }

            set_ghost_flags();
            set_boundary_flags();
        }

        /// <summary>
        /// Copy constructor
        /// </summary>
        public WindInflow(WindInflow old)
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
        }
    }
}
