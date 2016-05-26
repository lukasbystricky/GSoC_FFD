using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace FastFluidSolver
{
    public class PostProcessor
    {
        FluidSolver fs;
        Domain omega;

        /***************************************************************
         * Constructor
         **************************************************************/
        public PostProcessor(FluidSolver fs, Domain omega)
        {
            this.fs = fs;
            this.omega = omega;
        }

        /*****************************************************************************
        * Export data to a VTK file for visualization, based on file format guide here:
        * http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
         * 
         * Nx, Ny, Nz are the number of CELLS (not points) in each direction.
        ****************************************************************************/
        public void export_vtk(String fname, int Nx, int Ny, int Nz)
        {
            double hx = omega.length_x / Nx;
            double hy = omega.length_y / Ny;
            double hz = omega.length_z / Nz;

            double[, ,] p;
            double[, ,] u;
            double[, ,] v;
            double[, ,] w;

            interpolate_to_grid(Nx, Ny, Nz, out p, out u, out v, out w);

            using (StreamWriter sw = new StreamWriter(fname))
            {
                sw.WriteLine("# vtk DataFile Version 3.0");
                sw.WriteLine("Fast Fluid Dynamics data\n");
                sw.WriteLine("ASCII");
                sw.WriteLine("DATASET RECTILINEAR_GRID");
                sw.WriteLine("DIMENSIONS {0} {1} {2}", Nx + 1, Ny + 1, Nz + 1);
                sw.WriteLine("X_COORDINATES {0} double", (Nx + 1));
                for (int i = 0; i < Nx + 1; i++)
                {
                    sw.WriteLine("{0}", hx * i);
                }
                sw.WriteLine("Y_COORDINATES {0} double", (Ny + 1));
                for (int i = 0; i < Ny + 1; i++)
                {
                    sw.WriteLine("{0}", hy * i);
                }
                sw.WriteLine("Z_COORDINATES {0} double", (Nz + 1));
                for (int i = 0; i < Nz + 1; i++)
                {
                    sw.WriteLine("{0}", hz * i);
                }

                sw.WriteLine("POINT_DATA {0}", (Nx + 1) * (Ny + 1) * (Nz + 1));
                sw.WriteLine("VECTORS velocity double");
                for (int k = 0; k < Nz + 1; k++)
                {
                    for (int j = 0; j < Ny + 1; j++)
                    {
                        for (int i = 0; i < Nx + 1; i++)
                        {
                            sw.WriteLine("{0} {1} {2}", u[i, j, k], v[i, j, k], w[i, j, k]);
                        }
                    }
                }

                sw.WriteLine("SCALARS pressure double {0}", 1);
                sw.WriteLine("LOOKUP_TABLE default");
                for (int k = 0; k < Nz + 1; k++)
                {
                    for (int j = 0; j < Ny + 1; j++)
                    {
                        for (int i = 0; i < Nx + 1; i++)
                        {
                            sw.WriteLine("{0}", p[i, j, k]);
                        }
                    }
                }
            }
        }

        /************************************************************************************
         * Interpolates pressures and velocities from fs to a grid with Nx cells in x direction
         * Ny cells in y direction and Nz cells in z direction
         * 
         * Nx, Ny, Nz are the number of CELLS (not points) in each direction.
         * TO DO: implement in parallel 
         ***********************************************************************************/
        private void interpolate_to_grid(int Nx, int Ny, int Nz, out double[, ,] p_interp,
                        out double[, ,] u_interp, out double[, ,] v_interp, out double[, ,] w_interp)
        {
            p_interp = new double[Nx + 1, Ny + 1, Nz + 1];
            u_interp = new double[Nx + 1, Ny + 1, Nz + 1];
            v_interp = new double[Nx + 1, Ny + 1, Nz + 1];
            w_interp = new double[Nx + 1, Ny + 1, Nz + 1];

            double hx_interp = omega.length_x / Nx;
            double hy_interp = omega.length_y / Ny;
            double hz_interp = omega.length_z / Nz;

            double[] spacing_fs = new double[] { omega.hx, omega.hy, omega.hz };
            double[] coordinate = new double[3];

            for (int i = 0; i < Nx + 1; i++)
            {
                for (int j = 0; j < Ny + 1; j++)
                {
                    for (int k = 0; k < Nz + 1; k++)
                    {
                        coordinate[0] = i * hx_interp;
                        coordinate[1] = j * hy_interp;
                        coordinate[2] = k * hz_interp;

                        p_interp[i, j, k] = Utilities.trilinear_interpolation(coordinate, fs.p, 1, spacing_fs);
                        u_interp[i, j, k] = Utilities.trilinear_interpolation(coordinate, fs.u, 2, spacing_fs);
                        v_interp[i, j, k] = Utilities.trilinear_interpolation(coordinate, fs.v, 3, spacing_fs);
                        w_interp[i, j, k] = Utilities.trilinear_interpolation(coordinate, fs.w, 4, spacing_fs);
                    }
                }
            }
        }
    }
}

