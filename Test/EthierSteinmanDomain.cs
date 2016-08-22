using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/*
 * EthierSteinmanDomain.cs
 * Copyright 2016 Lukas Bystricky <lb13f@my.fsu.edu>
 *
 * This work is licensed under the GNU GPL license version 2 or later.
 */
 
namespace FastFluidSolver
{
    /// <summary>
    /// A cube domain with boundary conditions that give the exact solution described
    /// by Ethier and Steinman in the 1994 paper "Exact Fully 3D Navier-Stokes Solutions 
    /// for Benchmarking"
    /// </summary>
    public class EthierSteinmanDomain : Domain
    {
        private double nu, a, d;
        
        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="Nx">Number of cells in x direction (including ghost cells)</param>
        /// <param name="Ny">Number of cells in y direction (including ghost cells)</param>
        /// <param name="Nz">Number of cells in z direction (including ghost cells)</param>
        /// <param name="nu">fluid viscosity</param>
        /// <param name="a">coefficient a in boundary conditions/exact solutions<param>
        /// <param name="d">coefficient d in boundary conditions/exact solutions</param>
        public EthierSteinmanDomain(int Nx, int Ny, int Nz, double nu, double a, double d) 
        {
            this.Nx = Nx;
            this.Ny = Ny;
            this.Nz = Nz;

            this.nu = nu;
            this.a = a;
            this.d = d;

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
            set_ghost_flags();
        }

        /// <summary>
        /// Copy constructor
        /// </summary>
        public EthierSteinmanDomain(EthierSteinmanDomain old)
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

        /// <summary>
        /// Updates time dependent boundary conditions.
        /// </summary>
        /// <param name="t">time</param>
        public void update_boundary_conditions(double t)
        {
            /**************************************************************************************
            * u boundary conditions
            *************************************************************************************/
            for (int i = 0; i < boundary_u.GetLength(0); i++)
            {
                for (int j = 0; j < boundary_u.GetLength(1); j++)
                {
                    for (int k = 0; k < boundary_u.GetLength(2); k++)
                    {
                        double x = i*hx;
                        double y = (j - 0.5) * hy;
                        double z = (k - 0.5) * hz;

                        boundary_u[i, j, k] = -a * (Math.Exp(a * x) * Math.Sin(a * y + d * z) + 
                            Math.Exp(a * z) * Math.Cos(a * x + d * y)) * Math.Exp(-nu * Math.Pow(d, 2) * t);
                    }
                }
            }

            /**************************************************************************************
             * v boundary conditions
              *************************************************************************************/
            for (int i = 0; i < boundary_v.GetLength(0); i++)
            {
                for (int j = 0; j < boundary_v.GetLength(1); j++)
                {
                    for (int k = 0; k < boundary_v.GetLength(2); k++)
                    {
                        double x = (i - 0.5) * hx;
                        double y = j * hy;
                        double z = (k - 0.5) * hz;

                        boundary_v[i, j, k] = -a * (Math.Exp(a * y) * Math.Sin(a * z + d * x) +
                            Math.Exp(a * x) * Math.Cos(a * y + d * z)) * Math.Exp(-nu * Math.Pow(d, 2) * t);
                    }
                }
            }

            /**************************************************************************************
            * w boundary conditions
            *************************************************************************************/
            for (int i = 0; i < boundary_w.GetLength(0); i++)
            {
                for (int j = 0; j < boundary_w.GetLength(1); j++)
                {
                    for (int k = 0; k < boundary_w.GetLength(2); k++)
                    {
                        double x = (i - 0.5) * hx;
                        double y = (j - 0.5) * hy;
                        double z = k * hz;

                        boundary_w[i, j, k] = -a * (Math.Exp(a * z) * Math.Sin(a * x + d * y) +
                            Math.Exp(a * y) * Math.Cos(a * z + d * x)) * Math.Exp(-nu * Math.Pow(d, 2) * t);
                    }
                }
            }

        }

        /// <summary>
        /// Computes the exact solution at a point.
        /// </summary>
        /// <param name="coordinate">coorindate (x,y,z)</param>
        /// <param name="nu">fluid viscosity</param>
        /// <param name="t">time</param>
        /// <param name="u_exact">exact x component of velocity</param>
        /// <param name="v_exact">exact y component of velocity</param>
        /// <param name="w_exact">exact z component of velocity</param>
        /// <param name="p_exact">exact pressure</param>
        public override void exact_solution(double[] coordinate, double nu, double t, out double u_exact,
            out double v_exact, out double w_exact, out double p_exact)
        {
            double x = coordinate[0];
            double y = coordinate[1];
            double z = coordinate[2];

            u_exact = -a * (Math.Exp(a * x) * Math.Sin(a * y + d * z) +
                            Math.Exp(a * z) * Math.Cos(a * x + d * y)) * Math.Exp(-nu * Math.Pow(d, 2) * t);

            v_exact = -a * (Math.Exp(a * y) * Math.Sin(a * z + d * x) +
                            Math.Exp(a * x) * Math.Cos(a * y + d * z)) * Math.Exp(-nu * Math.Pow(d, 2) * t);

            w_exact = -a * (Math.Exp(a * z) * Math.Sin(a * x + d * y) +
                            Math.Exp(a * y) * Math.Cos(a * z + d * x)) * Math.Exp(-nu * Math.Pow(d, 2) * t);

            p_exact = (-Math.Pow(a, 2) / 2 * (Math.Exp(2 * a * x) + Math.Exp(2 * a * y) + Math.Exp(2 * a * z) +
                                    2 * Math.Sin(a * x + d * y) * Math.Cos(a * z + d * x) * Math.Exp(a * (y + z)) +
                                    2 * Math.Sin(a * y + d * z) * Math.Cos(a * x + d * y) * Math.Exp(a * (z + x)) +
                                    2 * Math.Sin(a * z + d * x) * Math.Cos(a * y + d * z) * Math.Exp(a * (x + y))) *
                                    Math.Exp(-2 * nu * Math.Pow(d, 2) * t));
        }
    }
}
