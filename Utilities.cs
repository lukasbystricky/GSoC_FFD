using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/*
 * Utilities.cs
 * Copyright 2016 Lukas Bystricky <lb13f@my.fsu.edu>
 *
 * This work is licensed under the GNU GPL license version 2 or later.
 */
 
namespace FastFluidSolver
{
    /// <summary>
    /// Ultilities that are common to multiple classes.
    /// </summary>
    public static class Utilities
    {
        const double EPS = 1e-12;

        /// <summary>
        /// Performs a trilinear interpolation. Method based on code in "Fluid Flow for 
        /// the Rest of Us: Tutorial of the Marker and Cell Method in Computer Graphics" 
        /// by Cline, Cardon and Egbert.
        /// </summary>
        /// <param name="x">Cell number in x direction</param>
        /// <param name="y">Cell number in y direction</param>
        /// <param name="z">Cell number in z direction</param>
        /// <param name="array">Array to interpolate</param>
        /// <returns>Interpolated value</returns>
        /// <remarks>Cell number here is the cell number on the underlying grid of the array. 
        /// They can be fractional and do not include ghost cells.</remarks>
        public static double trilinear_interpolation(double x, double y, double z,
                    double[, ,] array)
        {
            int imin = Math.Max(Math.Min((int)Math.Floor(x - EPS), array.GetLength(0) - 1), 0);
            int jmin = Math.Max(Math.Min((int)Math.Floor(y - EPS), array.GetLength(1) - 1), 0);
            int kmin = Math.Max(Math.Min((int)Math.Floor(z - EPS), array.GetLength(2) - 1), 0);

            int imax = Math.Max(Math.Min(imin + 1, array.GetLength(0) - 1), 0);
            int jmax = Math.Max(Math.Min(jmin + 1, array.GetLength(1) - 1), 0);
            int kmax = Math.Max(Math.Min(kmin + 1, array.GetLength(2) - 1), 0);

            return (imax - x) * (jmax - y) * (kmax - z) * array[imin, jmin, kmin] +
                    (x - imin) * (jmax - y) * (kmax - z) * array[imax, jmin, kmin] +
                    (imax - x) * (y - jmin) * (kmax - z) * array[imin, jmax, kmin] +
                    (x - imin) * (y - jmin) * (kmax - z) * array[imax, jmax, kmin] +
                    (imax - x) * (jmax - y) * (z - kmin) * array[imin, jmin, kmax] +
                    (x - imin) * (jmax - y) * (z - kmin) * array[imax, jmin, kmax] +
                    (imax - x) * (y - jmin) * (z - kmin) * array[imin, jmax, kmax] +
                    (x - imin) * (y - jmin) * (z - kmin) * array[imax, jmax, kmax];
        }

        /// <summary>
        /// Computes the pointwise L2 difference between 2 multidimensional arrays.
        /// </summary>
        /// <param name="x1">Array 1</param>
        /// <param name="x2">Array 2</param>
        /// <returns>Normalized pointwise L2 difference between array 1 and array 2.</returns>
        public static double compute_L2_difference(double[, ,] x1, double[, ,] x2)
        {
            double diff = 0;

            int Sx = x1.GetLength(0);
            int Sy = x1.GetLength(1);
            int Sz = x1.GetLength(2);

            for (int i = 0; i < Sx; i++)
            {
                for (int j = 0; j < Sy; j++)
                {
                    for (int k = 0; k < Sz; k++)
                    {
                        diff += Math.Pow(x1[i, j, k] - x2[i, j, k], 2);
                    }
                }
            }

            diff = Math.Sqrt(diff) / Math.Sqrt(Sx * Sy * Sz);
            return diff;
        }


        /// <summary>
        /// Checks if a point is inside the domain.
        /// </summary>
        /// <param name="coordinate">coordiate (x,y,z) of point</param>
        /// <param name="omega">Domain</param>
        /// <returns>True if point is inside fluid domain or on boundary, false if point is
        /// outside domain or inside and obstacle/ghost cell.</returns>
        public static bool in_domain(double[] coordinate, Domain omega)
        {
            double hx = omega.hx;
            double hy = omega.hy;
            double hz = omega.hz;

            //find domain cell that contains point
            int idomain_min = Math.Max(Math.Min((int)Math.Floor((1 - EPS) * (coordinate[0] + hx) / hx), omega.Nx - 1), 0);
            int jdomain_min = Math.Max(Math.Min((int)Math.Floor((1 - EPS) * (coordinate[1] + hy) / hy), omega.Ny - 1), 0);
            int kdomain_min = Math.Max(Math.Min((int)Math.Floor((1 - EPS) * (coordinate[2] + hz) / hz), omega.Nz - 1), 0);

            int idomain_max = Math.Max(Math.Min((int)Math.Floor((1 + EPS) * (coordinate[0] + hx) / hx), omega.Nx - 1), 0);
            int jdomain_max = Math.Max(Math.Min((int)Math.Floor((1 + EPS) * (coordinate[1] + hy) / hy), omega.Ny - 1), 0);
            int kdomain_max = Math.Max(Math.Min((int)Math.Floor((1 + EPS) * (coordinate[2] + hz) / hz), omega.Nz - 1), 0);


            List<int> possibleCelli = new List<int>();
            List<int> possibleCellj = new List<int>();
            List<int> possibleCellk = new List<int>();

            possibleCelli.Add(idomain_min);
            possibleCellj.Add(jdomain_min);
            possibleCellk.Add(kdomain_min);

            if (idomain_min != idomain_max)
            {
                possibleCelli.Add(idomain_max);
            }

            if (jdomain_min != jdomain_max)
            {
                possibleCellj.Add(jdomain_max);
            }

            if (kdomain_min != kdomain_max)
            {
                possibleCellk.Add(kdomain_max);
            }

            bool indomain = false;
            foreach (int i in possibleCelli)
            {
                foreach (int j in possibleCellj)
                {
                    foreach (int k in possibleCellk)
                    {
                        if (omega.obstacle_cells[i, j, k] == 0)
                        {
                            indomain = true;
                            break;
                        }
                    }
                }
            }

            return indomain;
        }

        /// <summary>
        /// Computes L2 and L-infinity error between FFD approximation and exact solution.
        /// </summary>
        /// <param name="fs">FluidSolver containing FFD solution at time t</param>
        /// <pram name="omega">Domain on which exact solution is given</pram>
        /// <param name="err_l2">Normalized pointwise L2 error</param>
        /// <param name="err_inf">Pointwise L-infinity error</param>
        /// <param name="t">Time</param>
        /// <param name="component">Component to evaluate</param>
        /// <remarks>The component to evaluate must be one of:
        ///     1 for pressure
        ///     2 for x component of velocity
        ///     3 for y component of velocity
        ///     4 for z component of velocity</remarks>
        public static void calculate_errors(FluidSolver fs, Domain omega, out double err_l2, out double err_inf,
            double t, int component)
        {
            DataExtractor de = new DataExtractor(omega, fs);

            double nu = fs.nu;

            int Nx = omega.Nx - 1;
            int Ny = omega.Ny - 1;
            int Nz = omega.Nz - 1;

            double[, ,] err_array = new double[Nx, Ny, Nz];
            double[, ,] comp_interp = new double[Nx, Ny, Nz];
            double[, ,] zeros = new double[Nx, Ny, Nz];

            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        double x = i * omega.hx;
                        double y = j * omega.hy;
                        double z = k * omega.hz;

                        double[] coordinate = new double[] { x, y, z };

                        double[] velocity_interp = de.get_velocity(x, y, z);

                        double u_exact, v_exact, w_exact, p_exact;

                        omega.exact_solution(coordinate, nu, t, out u_exact, out v_exact,
                                out w_exact, out p_exact);

                        switch (component)
                        {
                            case 1:
                                comp_interp[i, j, k] = de.get_pressure(x, y, z);
                                err_array[i, j, k] = Math.Abs(de.get_pressure(x, y, z) - p_exact);
                                break;

                            case 2:
                                comp_interp[i, j, k] = velocity_interp[0];
                                err_array[i, j, k] = Math.Abs(velocity_interp[0] - u_exact);
                                break;

                            case 3:
                                comp_interp[i, j, k] = velocity_interp[1];
                                err_array[i, j, k] = Math.Abs(velocity_interp[1] - v_exact);
                                break;

                            case 4:
                                comp_interp[i, j, k] = velocity_interp[2];
                                err_array[i, j, k] = Math.Abs(velocity_interp[2] - w_exact);
                                break;

                        }
                    }
                }
            }

            err_l2 = Utilities.compute_L2_difference(err_array, zeros); //L2 norm of errors
            err_inf = err_array.Cast<double>().Max();
        }
    }
}
