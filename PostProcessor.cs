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
        ****************************************************************************/
        public void export_vtk(String fname, int Nx, int Ny, int Nz)
        {

            double hx = omega.length_x / Nx;
            double hy = omega.length_y / Nx;
            double hz = omega.length_z / Nx;

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
                sw.WriteLine("DATASET STRUCTURED_GRID");
                sw.WriteLine("DIMENSIONS {0} {1} {2}", Nx, Ny, Nz);
                sw.WriteLine("POINTS {0} double", Nx * Ny * Nz);
                for (int i = 0; i < Nx; i++)
                {
                    for (int j = 0; j < Ny; j++)
                    {
                        for (int k = 0; k < Nz; k++)
                        {
                            sw.WriteLine("{0} {1} {2}", hx * i, hy * j, hz * k);
                        }
                    }
                }

                sw.WriteLine("POINT_DATA {0}", Nx * Ny * Nz);
                sw.WriteLine("VECTORS velocity double");
                for (int i = 0; i < Nx; i++)
                {
                    for (int j = 0; j < Ny; j++)
                    {
                        for (int k = 0; k < Nz; k++)
                        {
                            sw.WriteLine("{0} {1} {2}", fs.u[i, j, k], fs.v[i, j, k], fs.w[i, j, k]);
                        }
                    }
                }

                sw.WriteLine("SCALARS pressure double {0}", 1);
                sw.WriteLine("LOOKUP_TABLE default");
                for (int i = 0; i < Nz; i++)
                {
                    for (int j = 0; j < Ny; j++)
                    {
                        for (int k = 0; k < Nz; k++)
                        {
                            sw.WriteLine("{0}", fs.p[i, j, k]);
                        }
                    }
                }
            }
        }

        /************************************************************************************
         * Interpolates pressures and velocities from fs to a grid with Nx cells in x direction
         * Ny cells in y direction and Nz cells in z direction
         ***********************************************************************************/
        private void interpolate_to_grid(int Nx, int Ny, int Nz, out double[, ,] p_interp,
                        out double[, ,] u_interp, out double[, ,] v_interp, out double[, ,] w_interp)
        {
            p_interp = new double[Nx, Ny, Nz];
            u_interp = new double[Nx, Ny, Nz];
            v_interp = new double[Nx, Ny, Nz];
            w_interp = new double[Nx, Ny, Nz];

            double hx_interp = omega.length_x / Nx;
            double hy_interp = omega.length_y / Ny;
            double hz_interp = omega.length_z / Nz;

            int[] ncells_fs = new int[] { omega.Nx, omega.Ny, omega.Nz };
            double[] spacing_fs = new double[] { omega.hx, omega.hy, omega.hz };
            double[] coordinate = new double[3];

            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        coordinate[0] = i * hx_interp;
                        coordinate[1] = j * hy_interp;
                        coordinate[2] = k * hz_interp;

                        p_interp[i, j, k] = Utilities.trilinear_interpolation(coordinate, fs.p, 1, spacing_fs, ncells_fs);
                        u_interp[i, j, k] = Utilities.trilinear_interpolation(coordinate, fs.u, 2, spacing_fs, ncells_fs);
                        v_interp[i, j, k] = Utilities.trilinear_interpolation(coordinate, fs.v, 3, spacing_fs, ncells_fs);
                        w_interp[i, j, k] = Utilities.trilinear_interpolation(coordinate, fs.w, 4, spacing_fs, ncells_fs);
                    }
                }
            }
        }
    }
}

