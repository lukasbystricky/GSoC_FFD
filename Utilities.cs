using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolver
{    
    public static class Utilities
    {
        const double EPS = 1e-8;

        /*****************************************************************************
         * Performs a trilinear interpolation 
         * 
         * Method based on code in "Fluid Flow for the Rest of Us: Tutorial of the
         * Marker and Cell Method in Computer Graphics" by Cline, Cardon and Egbert
        ****************************************************************************/
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
