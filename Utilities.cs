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
        public static double trilinear_interpolation(double[] coordinate, double[, ,] array, int grid_type, 
            double[] spacings)
        {
            double hx = spacings[0];
            double hy = spacings[1];
            double hz = spacings[2];

            int Nx = array.GetLength(0);
            int Ny = array.GetLength(1);
            int Nz = array.GetLength(2);

            //Find grid cell that contains point
            int imin, jmin, kmin, imax, jmax, kmax;
            imin = imax = jmin = jmax = kmin = kmax = 0;

            double x, y, z;
            x = y = z = 0;

            switch (grid_type)
            {
                case 1://pressure grid

                    x = (coordinate[0]) % hx + 0.5;
                    y = (coordinate[1]) % hy + 0.5;
                    z = (coordinate[2]) % hz + 0.5;

                    imin = (int)Math.Floor(coordinate[0] / hx) + ((x > 0.5) ? 1 : 0);
                    jmin = (int)Math.Floor(coordinate[1] / hy) + ((y > 0.5) ? 1 : 0);
                    kmin = (int)Math.Floor(coordinate[2] / hz) + ((z > 0.5) ? 1 : 0);                    

                    break;

                case 2: //u velocity

                    x = (coordinate[0]) % hx;
                    y = (coordinate[1]) % hy + 0.5;
                    z = (coordinate[2]) % hz + 0.5;

                    imin = (int)Math.Floor(coordinate[0] / hx);
                    jmin = (int)Math.Floor(coordinate[1] / hy) + ((y > 0.5) ? 1 : 0);
                    kmin = (int)Math.Floor(coordinate[2] / hz) + ((z > 0.5) ? 1 : 0);

                    break;

                case 3: //v velocity

                    x = (coordinate[0]) % hx + 0.5;
                    y = (coordinate[1]) % hy;
                    z = (coordinate[2]) % hz + 0.5;

                    imin = (int)Math.Floor(coordinate[0] / hx) + ((x > 0.5) ? 1 : 0);
                    jmin = (int)Math.Floor(coordinate[1] / hy);
                    kmin = (int)Math.Floor(coordinate[2] / hz) + ((z > 0.5) ? 1 : 0);

                    break;

                case 4://w velocity

                    x = (coordinate[0]) % hx + 0.5;
                    y = (coordinate[1]) % hy + 0.5;
                    z = (coordinate[2]) % hz;

                    imin = (int)Math.Floor(coordinate[0] / hx) + ((x > 0.5) ? 1 : 0);
                    jmin = (int)Math.Floor(coordinate[1] / hy) + ((y > 0.5) ? 1 : 0);
                    kmin = (int)Math.Floor(coordinate[2] / hz);

                    break;
            }

            imax = Math.Min(imin + 1, Nx - 1);
            jmax = Math.Min(jmin + 1, Ny - 1);
            kmax = Math.Min(kmin + 1, Nz - 1);

            double[] corner_values = new double[8];            

            /*corner_values[0] = array[imin, jmin, kmin];
            corner_values[1] = array[imax, jmin, kmin];
            corner_values[2] = array[imin, jmin, kmax];
            corner_values[3] = array[imax, jmin, kmax];
            corner_values[4] = array[imin, jmax, kmin];
            corner_values[5] = array[imax, jmax, kmin];
            corner_values[6] = array[imin, jmax, kmax];
            corner_values[7] = array[imax, jmax, kmax];

            double c00 = corner_values[0] * (1 - x) + corner_values[1] * x;
            double c01 = corner_values[2] * (1 - x) + corner_values[3] * x;
            double c10 = corner_values[4] * (1 - x) + corner_values[5] * x;
            double c11 = corner_values[6] * (1 - x) + corner_values[7] * x;

            double c0 = c00 * (1 - y) + c10 * y;
            double c1 = c01 * (1 - y) + c11 * y;

            return c0 * (1 - z) + c1 * z;*/

            corner_values[0] = array[imin, jmin, kmin];
            corner_values[1] = array[imax, jmin, kmin];
            corner_values[2] = array[imin, jmax, kmin];
            corner_values[3] = array[imin, jmin, kmax];
            corner_values[4] = array[imax, jmin, kmax];
            corner_values[5] = array[imin, jmax, kmax];
            corner_values[6] = array[imax, jmax, kmin];
            corner_values[7] = array[imax, jmax, kmax];

            double v = corner_values[0] * (1 - x) * (1 - y) * (1 - z) + corner_values[1] * x * (1 - y) * (1 - z) +
                    corner_values[2] * (1 - x) * y * (1 - z) + corner_values[3] * (1 - x) * (1 - y) * z +
                    corner_values[4] * x * (1 - y) * z + corner_values[5] * (1 - x) * y * z +
                    corner_values[6] * x * y * (1 - z) + corner_values[7] * x * y * z;

            return v;
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

            diff = Math.Sqrt(diff) / Math.Sqrt(Sx * Sy * Sz);
            return diff;
        }
    }
}
