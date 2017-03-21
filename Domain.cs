using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/*
 * Domain.cs
 * Copyright 2016 Lukas Bystricky <lb13f@my.fsu.edu>
 *
 * This work is licensed under the GNU GPL license version 2 or later.
 */

namespace FastFluidSolver
{
    public class Domain
    {
        public int Nx, Ny, Nz;
        public double hx, hy, hz, length_x, length_y, length_z;

        public int[,,] boundary_cells { get; protected set; }   //flag to indicate if cell borders a boundary
        public int[,,] obstacle_cells { get; protected set; }   //flag to indicate if cell is part of an obstacle

        public int[,,] boundary_normal_x { get; protected set; } //flag to indicate if boundary at cell is normal to x direction
        public int[,,] boundary_normal_y { get; protected set; } //flag to indicate if boundary at cell is normal to y direction
        public int[,,] boundary_normal_z { get; protected set; } //flag to indicate if boundary at cell is normal to z direction

        public int[,,] outflow_boundary_x { get; protected set; }
        public int[,,] outflow_boundary_y { get; protected set; }
        public int[,,] outflow_boundary_z { get; protected set; }

        public double[,,] boundary_u { get; protected set; }   //x component of velocity at boundary
        public double[,,] boundary_v { get; protected set; }   //y component of velocity at boundary
        public double[,,] boundary_w { get; protected set; }   //z component of velocity at boundary

        public List<int[]> obstacle_list = new List<int[]>();
        public List<int[]> normal_x_list = new List<int[]>();
        public List<int[]> normal_y_list = new List<int[]>();
        public List<int[]> normal_z_list = new List<int[]>();

        /// <summary>
        /// Labels "ghost cells" outside domain as obstacles and flags cells adjacent to them as boundary cells.
        /// </summary>
        protected void set_ghost_flags()
        {
            //x = 0 and x = length_x boundaries
            for (int j = 0; j < Ny; j++)
            {
                for (int k = 0; k < Nz; k++)
                {
                    obstacle_cells[0, j, k] = 1;
                    obstacle_cells[Nx - 1, j, k] = 1;

                    obstacle_list.Add(new int[] { 0, j, k });
                    obstacle_list.Add(new int[] { Nx - 1, j, k });

                    boundary_cells[1, j, k] = 1;
                    boundary_cells[Nx - 2, j, k] = 1;

                    boundary_normal_x[1, j, k] = -1;
                    boundary_normal_x[Nx - 2, j, k] = 1;

                    normal_x_list.Add(new int[] { 1, j, k, -1 });
                    normal_x_list.Add(new int[] { Nx - 2, j, k, 1 });
                }
            }

            //y = 0 and y = length_y boundaries
            for (int i = 0; i < Nx; i++)
            {
                for (int k = 0; k < Nz; k++)
                {
                    obstacle_cells[i, 0, k] = 1;
                    obstacle_cells[i, Ny - 1, k] = 1;

                    obstacle_list.Add(new int[] { i, 0, k });
                    obstacle_list.Add(new int[] { i, Ny - 1, k });

                    boundary_cells[i, 1, k] = 1;
                    boundary_cells[i, Ny - 2, k] = 1;

                    boundary_normal_y[i, 1, k] = -1;
                    boundary_normal_y[i, Ny - 2, k] = 1;

                    normal_y_list.Add(new int[] { i, 1, k, -1 });
                    normal_y_list.Add(new int[] { i, Ny - 2, k, 1 });
                }
            }

            //z = 0 and z = length_z boundaries
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    obstacle_cells[i, j, 0] = 1;
                    obstacle_cells[i, j, Nz - 1] = 1;

                    obstacle_list.Add(new int[] { i, j, 0 });
                    obstacle_list.Add(new int[] { i, j, Nz - 1 });

                    boundary_cells[i, j, 1] = 1;
                    boundary_cells[i, j, Nz - 2] = 1;

                    boundary_normal_z[i, j, 1] = -1;
                    boundary_normal_z[i, j, Nz - 2] = 1;

                    normal_z_list.Add(new int[] { i, j, 1, -1 });
                    normal_z_list.Add(new int[] { i, j, Nz - 2, 1 });
                }
            }
        }

        /// <summary>
        /// Adds boundary flags to cells adjacent to obstacles.
        /// </summary>
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

                                normal_x_list.Add(new int[] { i + 1, j, k, -1 });
                            }

                            if (obstacle_cells[i - 1, j, k] == 0)
                            {
                                boundary_cells[i - 1, j, k] = 1;
                                boundary_normal_x[i - 1, j, k] = 1;

                                normal_x_list.Add(new int[] { i - 1, j, k, 1 });
                            }

                            if (obstacle_cells[i, j + 1, k] == 0)
                            {
                                boundary_cells[i, j + 1, k] = 1;
                                boundary_normal_y[i, j + 1, k] = -1;

                                normal_y_list.Add(new int[] { i, j + 1, k, -1 });
                            }

                            if (obstacle_cells[i, j - 1, k] == 0)
                            {
                                boundary_cells[i, j - 1, k] = 1;
                                boundary_normal_y[i, j - 1, k] = 1;

                                normal_y_list.Add(new int[] { i, j - 1, k, 1 });
                            }

                            if (obstacle_cells[i, j, k + 1] == 0)
                            {
                                boundary_cells[i, j, k + 1] = 1;
                                boundary_normal_z[i, j, k + 1] = -1;

                                normal_z_list.Add(new int[] { i, j, k + 1, -1 });
                            }

                            if (obstacle_cells[i, j, k - 1] == 0)
                            {
                                boundary_cells[i, j, k - 1] = 1;
                                boundary_normal_z[i, j, k - 1] = 1;

                                normal_z_list.Add(new int[] { i, j, k - 1, 1 });
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Adds an obstacle to the domain.
        /// </summary>
        /// <param name="xmin">minimum x coordinate of obstacle</param>
        /// <param name="xmax">maximum x coordinate of obstacle</param>
        /// <param name="ymin">minimum y coordinate of obstacle</param>
        /// <param name="ymax">maximum y coordinate of obstacle</param>
        /// <param name="zmin">minimum z coordinate of obstacle</param>
        /// <param name="zmax">maximum z coordinate of obstacle</param>
        public void add_obstacle(double xmin, double xmax, double ymin, double ymax,
                double zmin, double zmax)
        {
            int i_start = (int)Math.Floor(xmin * (Nx - 2) / length_x);
            int i_end = (int)Math.Floor(xmax * (Nx - 2) / length_x);
            int j_start = (int)Math.Floor(ymin * (Ny - 2) / length_y);
            int j_end = (int)Math.Floor(ymax * (Ny - 2) / length_y);
            int k_start = (int)Math.Floor(zmin * (Nz - 2) / length_z);
            int k_end = (int)Math.Floor(zmax * (Nz - 2) / length_z);

            for (int i = i_start; i < i_end; i++)
            {
                for (int j = j_start; j < j_end; j++)
                {
                    for (int k = k_start; k < k_end; k++)
                    {
                        obstacle_cells[i + 1, j + 1, k + 1] = 1;
                        //obstacle_list.Add(new int[] { i, j, k });
                    }
                }
            }

            set_boundary_flags();
        }

        /// <summary>
        /// Computes the exact solution at a point if known, the default returns 0 for everything.
        /// </summary>
        /// <param name="coordinate">coorindate (x,y,z)</param>
        /// <param name="nu">fluid viscosity</param>
        /// <param name="t">time</param>
        /// <param name="u_exact">exact x component of velocity</param>
        /// <param name="v_exact">exact y component of velocity</param>
        /// <param name="w_exact">exact z component of velocity</param>
        /// <param name="p_exact">exact pressure</param>
        virtual public void exact_solution(double[] coordinate, double nu, double t, out double u_exact,
            out double v_exact, out double w_exact, out double p_exact)
        {
            u_exact = 0;
            v_exact = 0;
            w_exact = 0;
            p_exact = 0;
        }
    }
}