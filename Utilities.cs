using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolver
{    
    public static class Utilities
    {
        const double EPS = 1e-6;

        /*****************************************************************************
         * Performs a trilinear interpolation 
        ****************************************************************************/
        public static double trilinear_interpolation(double[] coordinate, double[, ,] array, int grid_type, Domain omega)
        {
            double hx = omega.hx;
            double hy = omega.hy;
            double hz = omega.hz;

            int Nx = array.GetLength(0);
            int Ny = array.GetLength(1);
            int Nz = array.GetLength(2);

            //check if point inside obstacle
            if (!in_domain(coordinate, omega))
            {
                return 0;
            }

            //Find grid cell that contains point
            int imin, jmin, kmin, imax, jmax, kmax;
            imin = imax = jmin = jmax = kmin = kmax = 0;

            double x, y, z;
            x = y = z = 0;

            switch (grid_type)
            {
                case 1://cell centred grid

                    x = coordinate[0] - Math.Floor(coordinate[0] / hx) * hx + 0.5;
                    y = coordinate[1] - Math.Floor(coordinate[1] / hy) * hz + 0.5;
                    z = coordinate[2] - Math.Floor(coordinate[2] / hz) * hz + 0.5;

                    imin = (int)Math.Floor(coordinate[0] / hx) + ((x > 0.5) ? 1 : 0);
                    jmin = (int)Math.Floor(coordinate[1] / hy) + ((y > 0.5) ? 1 : 0);
                    kmin = (int)Math.Floor(coordinate[2] / hz) + ((z > 0.5) ? 1 : 0);                    

                    break;

                case 2: //u velocity

                    x = coordinate[0] - Math.Floor(coordinate[0] / hx) * hx;
                    y = coordinate[1] - Math.Floor(coordinate[1] / hy) * hz + 0.5;
                    z = coordinate[2] - Math.Floor(coordinate[2] / hz) * hz + 0.5;

                    imin = (int)Math.Floor(coordinate[0] / hx);
                    jmin = (int)Math.Floor(coordinate[1] / hy) + ((y > 0.5) ? 1 : 0);
                    kmin = (int)Math.Floor(coordinate[2] / hz) + ((z > 0.5) ? 1 : 0);

                    break;

                case 3: //v velocity

                    x = coordinate[0] - Math.Floor(coordinate[0] / hx) * hx + 0.5;
                    y = coordinate[1] - Math.Floor(coordinate[1] / hy) * hz;
                    z = coordinate[2] - Math.Floor(coordinate[2] / hz) * hz + 0.5;

                    imin = (int)Math.Floor(coordinate[0] / hx) + ((x > 0.5) ? 1 : 0);
                    jmin = (int)Math.Floor(coordinate[1] / hy);
                    kmin = (int)Math.Floor(coordinate[2] / hz) + ((z > 0.5) ? 1 : 0);

                    break;

                case 4://w velocity

                    x = coordinate[0] - Math.Floor(coordinate[0] / hx) * hx + 0.5;
                    y = coordinate[1] - Math.Floor(coordinate[1] / hy) * hz + 0.5;
                    z = coordinate[2] - Math.Floor(coordinate[2] / hz) * hz;

                    imin = (int)Math.Floor(coordinate[0] / hx) + ((x > 0.5) ? 1 : 0);
                    jmin = (int)Math.Floor(coordinate[1] / hy) + ((y > 0.5) ? 1 : 0);
                    kmin = (int)Math.Floor(coordinate[2] / hz);

                    break;
            }

            imax = Math.Min(imin + 1, Nx - 1);
            jmax = Math.Min(jmin + 1, Ny - 1);
            kmax = Math.Min(kmin + 1, Nz - 1);

            double[] corner_values = new double[8];            

            corner_values[0] = array[imin, jmin, kmin];
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

            return c0 * (1 - z) + c1 * z;

           /* corner_values[0] = array[imin, jmin, kmin];
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

            return v;*/
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

        /***********************************************************************
         * Checks if a point is inside the domain omega. Returns true if point
         * is inside domain or on a boundary, false otherwise
         ***********************************************************************/
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
            int jdomain_max = Math.Max(Math.Min((int)Math.Floor((1 + EPS) * (coordinate[1] + hy) / hy), omega.Nx - 1), 0);
            int kdomain_max = Math.Max(Math.Min((int)Math.Floor((1 + EPS) * (coordinate[2] + hz) / hz), omega.Nx - 1), 0);


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
    }
}
