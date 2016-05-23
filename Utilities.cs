using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolver
{
    public static class Utilities
    {
        /*****************************************************************************
         * Performs a trilinear interpolation 
        ****************************************************************************/
        public static double trilinear_interpolation(double[] coordinate, double[, ,] array, int grid_type, double[] spacing)
        {
            double hx = spacing[0];
            double hy = spacing[1];
            double hz = spacing[2];

            int i = (int)Math.Floor(coordinate[0] / hx);
            int j = (int)Math.Floor(coordinate[1] / hy);
            int k = (int)Math.Floor(coordinate[2] / hz);

            int imin, imax, jmin, jmax, kmin, kmax;
            double[] corner_values = new double[8];

            // interpolate pressure values
            // pressure values are given at centre of cell
            double[] centre = find_centre(i, j, k, grid_type);
            find_bounding_indices(centre[0], coordinate[0], i, out imin, out imax);
            find_bounding_indices(centre[1], coordinate[1], j, out jmin, out jmax);
            find_bounding_indices(centre[k], coordinate[k], k, out kmin, out kmax);

            corner_values[0] = array[imin, jmin, kmin];
            corner_values[1] = array[imax, jmin, kmin];
            corner_values[2] = array[imin, jmin, kmax];
            corner_values[3] = array[imax, jmin, kmax];
            corner_values[4] = array[imin, jmax, kmin];
            corner_values[5] = array[imax, jmax, kmin];
            corner_values[6] = array[imin, jmax, kmax];
            corner_values[7] = array[imax, jmax, kmax];

            double x = coordinate[0] - centre[0];
            double y = coordinate[1] - centre[1];
            double z = coordinate[2] - centre[2];

            double c00 = corner_values[0] * (hx - x) + corner_values[1] * x;
            double c01 = corner_values[2] * (hx - x) + corner_values[3] * x;
            double c10 = corner_values[4] * (hx - x) + corner_values[5] * x;
            double c11 = corner_values[6] * (hx - x) + corner_values[7] * x;

            double c0 = c00 * (hy - y) + c10 * y;
            double c1 = c01 * (hy - y) + c11 * y;

            return c0 * (hz - z) + c1 * z;
        }

        /**************************************************************************
         * Computes the normalized L2 difference between two 3 dimensional arrays
         **************************************************************************/
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

            diff = Math.Sqrt(diff) / (Sx * Sy * Sz);
            return diff;
        }

        /***************************************************************************
         * Find cartestian coordinates of centre of cell i, j, k
         **************************************************************************/
        public static double[] find_centre(int i, int j, int k, int grid_type)
        {
            double[] centre = new double[3];

            switch (grid_type)
            {
                case 1:
                    centre[0] = Math.Floor(i / hx) + hx / 2;
                    centre[1] = Math.Floor(j / hy) + hy / 2;
                    centre[2] = Math.Floor(k / hz) + hz / 2;
                    break;

                case 2:
                    centre[0] = Math.Floor(i / hx);
                    centre[1] = Math.Floor(j / hy) + hy / 2;
                    centre[2] = Math.Floor(k / hz) + hz / 2;
                    break;

                case 3:
                    centre[0] = Math.Floor(i / hx) + hx / 2;
                    centre[1] = Math.Floor(j / hy);
                    centre[2] = Math.Floor(k / hz) + hz / 2;
                    break;

                case 4:
                    centre[0] = Math.Floor(i / hx) + hx / 2;
                    centre[1] = Math.Floor(j / hy) + hy / 2;
                    centre[2] = Math.Floor(k / hz);
                    break;

            }

            return centre;
        }

        /***********************************************************************
         * Find the bounding indices for a linear interpolation
         **********************************************************************/
        public static void find_bounding_indices(double centre, double coordinate, int i, out int imin, out int imax)
        {
            if (centre > coordinate)
            {
                imin = i - 1;
                imax = i;
            }
            else
            {
                imin = i;
                imax = i + 1;
            }
        }
    }
}
