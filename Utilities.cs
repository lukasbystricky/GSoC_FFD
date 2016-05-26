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

            double[] grid_offset = new double[3];
            switch (grid_type)
            {
                case 1://pressure grid
                    grid_offset[0] = -0.5;
                    grid_offset[1] = -0.5;
                    grid_offset[2] = -0.5;

                    imin = Math.Min((int)Math.Floor(coordinate[0] / hx), Nx - 2);
                    jmin = Math.Min((int)Math.Floor(coordinate[1] / hy), Ny - 2);
                    kmin = Math.Min((int)Math.Floor(coordinate[2] / hz), Nz - 2);

                    break;

                case 2: //u velocity
                    grid_offset[1] = -0.5;
                    grid_offset[2] = -0.5;

                    imin = Math.Min((int)Math.Floor(coordinate[0] / hx), Nx - 2);
                    jmin = Math.Min((int)Math.Floor((coordinate[1] + 0.5 * hy) / hy), Ny - 2);
                    kmin = Math.Min((int)Math.Floor((coordinate[2] + 0.5 * hz) / hz), Nz - 2);

                    break;

                case 3: //v velocity
                    grid_offset[0] = -0.5;
                    grid_offset[2] = -0.5;

                    imin = Math.Min((int)Math.Floor((coordinate[0] + 0.5 * hx) / hx), Nx - 2);
                    jmin = Math.Min((int)Math.Floor(coordinate[1] / hy), Ny - 2);
                    kmin = Math.Min((int)Math.Floor((coordinate[2] + 0.5 * hz) / hz), Nz - 2);

                    break;

                case 4://w velocity
                    grid_offset[0] = -0.5;
                    grid_offset[1] = -0.5;

                    imin = Math.Min((int)Math.Floor((coordinate[0] + 0.5 * hx)/ hx), Nx - 2);
                    jmin = Math.Min((int)Math.Floor((coordinate[1] + 0.5 * hy)/ hy), Ny - 2);
                    kmin = Math.Min((int)Math.Floor(coordinate[2] / hz), Nz - 2);

                    break;
            }
            
            double[] corner_values = new double[8];

            imax = imin + 1;
            jmax = jmin + 1;
            kmax = kmin + 1;

            corner_values[0] = array[imin, jmin, kmin];
            corner_values[1] = array[imax, jmin, kmin];
            corner_values[2] = array[imin, jmin, kmax];
            corner_values[3] = array[imax, jmin, kmax];
            corner_values[4] = array[imin, jmax, kmin];
            corner_values[5] = array[imax, jmax, kmin];
            corner_values[6] = array[imin, jmax, kmax];
            corner_values[7] = array[imax, jmax, kmax];           

            double x = (coordinate[0] - hx * imin) / hx - grid_offset[0];
            double y = (coordinate[1] - hy * jmin) / hy - grid_offset[1];
            double z = (coordinate[2] - hz * kmin) / hz - grid_offset[2];

            double c00 = corner_values[0] * (1 - x) + corner_values[1] * x;
            double c01 = corner_values[2] * (1 - x) + corner_values[3] * x;
            double c10 = corner_values[4] * (1 - x) + corner_values[5] * x;
            double c11 = corner_values[6] * (1 - x) + corner_values[7] * x;

            double c0 = c00 * (1 - y) + c10 * y;
            double c1 = c01 * (1 - y) + c11 * y;

            return c0 * (1 - z) + c1 * z;
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
