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
        public void export_vtk(String fname, int Nx, int Ny, int Nz, double time)
        {
            int Nx_points = Nx + 1;
            int Ny_points = Ny + 1;
            int Nz_points = Nz + 1;

            double hx = omega.length_x / Nx;
            double hy = omega.length_y / Ny;
            double hz = omega.length_z / Nz;

            double[, ,] p_interp;
            double[, ,] u_interp;
            double[, ,] v_interp;
            double[, ,] w_interp;

            interpolate_to_grid(Nx, Ny, Nz, out p_interp, out u_interp, out v_interp, out w_interp);

            using (StreamWriter sw = new StreamWriter(fname))
            {
                sw.WriteLine("# vtk DataFile Version 3.0");
                sw.WriteLine("Fast Fluid Dynamics data\n");
                sw.WriteLine("ASCII");
                sw.WriteLine("DATASET RECTILINEAR_GRID");
                sw.WriteLine("FIELD FieldData 1");
                sw.WriteLine("TIME 1 1 double");
                sw.WriteLine("{0}", time);
                sw.WriteLine("DIMENSIONS {0} {1} {2}", Nx_points, Ny_points, Nz_points);
                sw.WriteLine("X_COORDINATES {0} double", Nx_points);

                for (int i = 0; i < Nx_points; i++)
                {
                    sw.WriteLine("{0}", hx * i);
                }
                sw.WriteLine("Y_COORDINATES {0} double", Ny_points);
                for (int i = 0; i < Ny_points; i++)
                {
                    sw.WriteLine("{0}", hy * i);
                }
                sw.WriteLine("Z_COORDINATES {0} double", Nz_points);
                for (int i = 0; i < Nz_points; i++)
                {
                    sw.WriteLine("{0}", hz * i);
                }

               sw.WriteLine("POINT_DATA {0}", Nx_points * Ny_points * Nz_points);
               sw.WriteLine("VECTORS velocity double");
               for (int k = 0; k < Nz_points; k++)
               {
                   for (int j = 0; j < Ny_points; j++)
                   {
                       for (int i = 0; i < Nx_points; i++)
                       {
                           sw.WriteLine("{0} {1} {2}", u_interp[i, j, k], v_interp[i, j, k], w_interp[i, j, k]);
                       }
                   }
               }
                
                sw.WriteLine("SCALARS pressure double {0}", 1);
                sw.WriteLine("LOOKUP_TABLE default");
                for (int k = 0; k < Nz_points; k++)
                {
                    for (int j = 0; j < Ny_points; j++)
                    {
                        for (int i = 0; i < Nx_points; i++)
                        {
                            sw.WriteLine("{0}", p_interp[i, j, k]);
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

            double[] coordinate = new double[3];

            for (int i = 0; i < Nx + 1; i++)
            {
                for (int j = 0; j < Ny + 1; j++)
                {
                    for (int k = 0; k < Nz + 1; k++)
                    {
                        u_interp[i, j, k] = Utilities.trilinear_interpolation(i, j + 0.5, k + 0.5, fs.u);
                        v_interp[i, j, k] = Utilities.trilinear_interpolation(i + 0.5, j, k + 0.5, fs.v);
                        w_interp[i, j, k] = Utilities.trilinear_interpolation(i + 0.5, j + 0.5, k, fs.w);
                        p_interp[i, j, k] = Utilities.trilinear_interpolation(i, j, k, fs.p);
                    }
                }
            }
        }
    }
}

