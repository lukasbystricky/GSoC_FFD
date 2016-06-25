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
        public void export_data_vtk(String fname, int Nx, int Ny, int Nz, double time)
        {
            double hx = omega.length_x / Nx;
            double hy = omega.length_y / Ny;
            double hz = omega.length_z / Nz;

            double[, ,] u_interp;
            double[, ,] v_interp;
            double[, ,] w_interp;

            interpolate_to_grid(Nx, Ny, Nz, out u_interp, out v_interp, out w_interp);

            using (StreamWriter sw = new StreamWriter(fname))
            {
                sw.WriteLine("# vtk DataFile Version 3.0");
                sw.WriteLine("Fast Fluid Dynamics data\n");
                sw.WriteLine("ASCII");
                sw.WriteLine("DATASET RECTILINEAR_GRID");
                sw.WriteLine("FIELD FieldData 1");
                sw.WriteLine("TIME 1 1 double");
                sw.WriteLine("{0}", time);
                sw.WriteLine("DIMENSIONS {0} {1} {2}", Nx + 1, Ny + 1, Nz + 1);
                sw.WriteLine("X_COORDINATES {0} double", Nx + 1);

                for (int i = 0; i < Nx + 1; i++)
                {
                    sw.WriteLine("{0}", hx * i);
                }
                sw.WriteLine("Y_COORDINATES {0} double", Ny + 1);
                for (int i = 0; i < Ny + 1; i++)
                {
                    sw.WriteLine("{0}", hy * i);
                }
                sw.WriteLine("Z_COORDINATES {0} double", Nz + 1);
                for (int i = 0; i < Nz + 1; i++)
                {
                    sw.WriteLine("{0}", hz * i);
                }

               sw.WriteLine("CELL_DATA {0}", Nx * Ny * Nz);
               sw.WriteLine("VECTORS velocity double");
               for (int k = 0; k < Nz; k++)
               {
                   for (int j = 0; j < Ny; j++)
                   {
                       for (int i = 0; i < Nx; i++)
                       {
                           sw.WriteLine("{0} {1} {2}", u_interp[i, j, k], v_interp[i, j, k], w_interp[i, j, k]);
                       }
                   }
               }
                
                sw.WriteLine("SCALARS pressure double {0}", 1);
                sw.WriteLine("LOOKUP_TABLE default");
                for (int k = 1; k < Nz + 1; k++)
                {
                    for (int j = 1; j < Ny + 1; j++)
                    {
                        for (int i = 1; i < Nx + 1; i++)
                        {
                            sw.WriteLine("{0}", fs.p[i, j, k]);
                        }
                    }
                }
            }
        }

        /*****************************************************************************
            * Export geometry to a VTK file for visualization
        *****************************************************************************/
       public void export_geometry_vtk(String fname, double time)
       {
           double hx = omega.length_x / (omega.Nx - 2);
           double hy = omega.length_y / (omega.Ny - 2);
           double hz = omega.length_z / (omega.Nz - 2);

           using (StreamWriter sw = new StreamWriter(fname))
           {
               sw.WriteLine("# vtk DataFile Version 3.0");
               sw.WriteLine("Fast Fluid Dynamics geometry\n");
               sw.WriteLine("ASCII");
               sw.WriteLine("DATASET RECTILINEAR_GRID");
               sw.WriteLine("FIELD FieldData 1");
               sw.WriteLine("TIME 1 1 double");
               sw.WriteLine("{0}", time);
               sw.WriteLine("DIMENSIONS {0} {1} {2}", omega.Nx + 1, omega.Ny + 1, omega.Nz + 1);
               sw.WriteLine("X_COORDINATES {0} double", omega.Nx + 1);

               for (int i = 0; i < omega.Nx + 1; i++)
               {
                   sw.WriteLine("{0}", hx * (i - 1));
               }
               sw.WriteLine("Y_COORDINATES {0} double", omega.Ny + 1);
               for (int i = 0; i < omega.Ny + 1; i++)
               {
                   sw.WriteLine("{0}", hy * (i - 1));
               }
               sw.WriteLine("Z_COORDINATES {0} double", omega.Nz + 1);
               for (int i = 0; i < omega.Nz + 1; i++)
               {
                   sw.WriteLine("{0}", hz * (i - 1));
               }

               sw.WriteLine("CELL_DATA {0}", omega.Nx * omega.Ny * omega.Nz);
               sw.WriteLine("SCALARS boundary_cells double {0}", 1);
               sw.WriteLine("LOOKUP_TABLE default");
               for (int k = 0; k < omega.Nz; k++)
               {
                   for (int j = 0; j < omega.Ny; j++)
                   {
                       for (int i = 0; i < omega.Nx; i++)
                       {
                           sw.Write("{0} ", (double) omega.boundary_cells[i, j, k]);
                       }
                   }
               }
               sw.WriteLine("SCALARS obstacle_cells double {0}", 1);
               sw.WriteLine("LOOKUP_TABLE default");
               for (int k = 0; k < omega.Nz; k++)
               {
                   for (int j = 0; j < omega.Ny; j++)
                   {
                       for (int i = 0; i < omega.Nx; i++)
                       {
                           sw.Write("{0} ", (double) omega.obstacle_cells[i, j, k]);
                       }
                   }
               }

               sw.WriteLine("SCALARS boundary_normal_x double {0}", 1);
               sw.WriteLine("LOOKUP_TABLE default");
               for (int k = 0; k < omega.Nz; k++)
               {
                   for (int j = 0; j < omega.Ny; j++)
                   {
                       for (int i = 0; i < omega.Nx; i++)
                       {
                           sw.Write("{0} ", (double) omega.boundary_normal_x[i, j, k]);
                       }
                   }
               }

               sw.WriteLine("SCALARS boundary_normal_y double {0}", 1);
               sw.WriteLine("LOOKUP_TABLE default");
               for (int k = 0; k < omega.Nz; k++)
               {
                   for (int j = 0; j < omega.Ny; j++)
                   {
                       for (int i = 0; i < omega.Nx; i++)
                       {
                           sw.Write("{0} ", (double) omega.boundary_normal_y[i, j, k]);
                       }
                   }
               }

               sw.WriteLine("SCALARS boundary_normal_z double {0}", 1);
               sw.WriteLine("LOOKUP_TABLE default");
               for (int k = 0; k < omega.Nz; k++)
               {
                   for (int j = 0; j < omega.Ny; j++)
                   {
                       for (int i = 0; i < omega.Nx; i++)
                       {
                           sw.Write("{0} ", (double) omega.boundary_normal_z[i, j, k]);
                       }
                   }
               }

               sw.WriteLine("SCALARS outflow_cells double {0}", 1);
               sw.WriteLine("LOOKUP_TABLE default");
               for (int k = 0; k < omega.Nz; k++)
               {
                   for (int j = 0; j < omega.Ny; j++)
                   {
                       for (int i = 0; i < omega.Nx; i++)
                       {
                           sw.Write("{0} ", (double) omega.outflow_cells[i, j, k]);
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
        private void interpolate_to_grid(int Nx, int Ny, int Nz,
                        out double[, ,] u_interp, out double[, ,] v_interp, out double[, ,] w_interp)
        {
            u_interp = new double[Nx, Ny, Nz];
            v_interp = new double[Nx, Ny, Nz];
            w_interp = new double[Nx, Ny, Nz];

            double hx_interp = omega.length_x / Nx;
            double hy_interp = omega.length_y / Ny;
            double hz_interp = omega.length_z / Nz;

            double[] coordinate = new double[3];

            for (int i = 1; i < Nx + 1; i++)
            {
                for (int j = 1; j < Ny + 1; j++)
                {
                    for (int k = 1; k < Nz + 1; k++)
                    {
                        if (omega.obstacle_cells[i,j,k] == 0)
                        {
                            u_interp[i - 1, j - 1, k - 1] = Utilities.trilinear_interpolation(i + 1, j, k, fs.u);
                            v_interp[i - 1, j - 1, k - 1] = Utilities.trilinear_interpolation(i, j - 0.5, k, fs.v);
                            w_interp[i - 1, j - 1, k - 1] = Utilities.trilinear_interpolation(i, j, k - 0.5, fs.w);                          
                        }
                    }
                }
            }
        }
    }
}

