using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolver
{
    public class Domain
    {
        public int Nx, Ny, Nz;
        public double hx, hy, hz, length_x, length_y, length_z;

        public int[, ,] boundary_cells { get; protected set; }   //flag to indicate if cell borders a boundary
        public int[, ,] obstacle_cells { get; protected set; }   //flag to indicate if cell is part of an obstacle

        public int[, ,] boundary_normal_x { get; protected set; } //flag to indicate if boundary at cell is normal to x direction
        public int[, ,] boundary_normal_y { get; protected set; } //flag to indicate if boundary at cell is normal to y direction
        public int[, ,] boundary_normal_z { get; protected set; } //flag to indicate if boundary at cell is normal to z direction

        public int[, ,] outflow_boundary_x { get; protected set; }
        public int[, ,] outflow_boundary_y { get; protected set; }
        public int[, ,] outflow_boundary_z { get; protected set; }

        public double[, ,] boundary_u { get; protected set; }   //x component of velocity at boundary
        public double[, ,] boundary_v { get; protected set; }   //y component of velocity at boundary
        public double[, ,] boundary_w { get; protected set; }   //z component of velocity at boundary  

        /*********************************************************************************************
         * add flags indicating obstacle (ghost) cells around entire domain
         ********************************************************************************************/
        protected void set_ghost_flags()
        {
            //x = 0 and x = length_x boundaries
            for (int j = 0; j < Ny; j++)
            {
                for (int k = 0; k < Nz; k++)
                {
                    obstacle_cells[0, j, k] = 1;
                    obstacle_cells[Nx - 1, j, k] = 1;

                    boundary_cells[1, j, k] = 1;
                    boundary_cells[Nx - 2, j, k] = 1;

                    boundary_normal_x[1, j, k] = -1;
                    boundary_normal_x[Nx - 2, j, k] = 1;
                }
            }

            //y = 0 and y = length_y boundaries
            for (int i = 0; i < Nx; i++)
            {
                for (int k = 0; k < Nz; k++)
                {
                    obstacle_cells[i, 0, k] = 1;
                    obstacle_cells[i, Ny - 1, k] = 1;

                    boundary_cells[i, 1, k] = 1;
                    boundary_cells[i, Ny - 2, k] = 1;

                    boundary_normal_y[i, 1, k] = -1;
                    boundary_normal_y[i, Ny - 2, k] = 1;
                }
            }

            //z = 0 and z = length_z boundaries
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    obstacle_cells[i, j, 0] = 1;
                    obstacle_cells[i, j, Nz - 1] = 1;

                    boundary_cells[i, j, 1] = 1;
                    boundary_cells[i, j, Nz - 2] = 1;

                    boundary_normal_z[i, j, 1] = -1;
                    boundary_normal_z[i, j, Nz - 2] = 1;
                }
            }
        }

        /*******************************************************************************************
         * add boundary flags and normal flags on cells bordering obstacles
         ******************************************************************************************/
        protected void set_boundary_flags()
        {
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    for (int k = 1; k < Nz - 1; k++)
                    {
                        if (obstacle_cells[i, j, k] == 1)
                        {
                            if (obstacle_cells[i + 1, j, k] == 0)
                            {
                                boundary_cells[i + 1, j, k] = 1;
                                boundary_normal_x[i + 1, j, k] = -1;
                            }

                            if (obstacle_cells[i - 1, j, k] == 0)
                            {
                                boundary_cells[i - 1, j, k] = 1;
                                boundary_normal_x[i - 1, j, k] = 1;
                            }

                            if (obstacle_cells[i, j + 1, k] == 0)
                            {
                                boundary_cells[i, j + 1, k] = 1;
                                boundary_normal_y[i, j + 1, k] = -1;
                            }

                            if (obstacle_cells[i, j - 1, k] == 0)
                            {
                                boundary_cells[i, j - 1, k] = 1;
                                boundary_normal_y[i, j - 1, k] = 1;
                            }

                            if (obstacle_cells[i, j, k + 1] == 0)
                            {
                                boundary_cells[i, j, k + 1] = 1;
                                boundary_normal_z[i, j, k + 1] = -1;
                            }

                            if (obstacle_cells[i, j, k - 1] == 0)
                            {
                                boundary_cells[i, j, k - 1] = 1;
                                boundary_normal_z[i, j, k - 1] = 1;
                            }
                        }
                    }
                }
            }
        }
    }
}
