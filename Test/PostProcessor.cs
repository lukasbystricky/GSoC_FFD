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
        DataExtractor de;

        /***************************************************************
         * Constructor
         **************************************************************/
        public PostProcessor(FluidSolver fs, Domain omega)
        {
            this.fs = fs;
            this.omega = omega;

            de = new DataExtractor(omega, fs);
        }

        /*****************************************************************************
        * Export data to a VTK file for visualization, based on file format guide here:
        * http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
        * 
        * Nx, Ny, Nz are the number of CELLS (not points) in each direction.
        ****************************************************************************/
        public void export_data_vtk(String fname, double time, bool cell_centred)
        {

            int Nx, Ny, Nz;

            double hx = omega.hx;
            double hy = omega.hy;
            double hz = omega.hz;

            double[, ,] u_interp;
            double[, ,] v_interp;
            double[, ,] w_interp;
            double[, ,] p_interp;
            double[, ,] div_interp;

            if (cell_centred)
            {
                Nx = omega.Nx - 2;
                Ny = omega.Ny - 2;
                Nz = omega.Nz - 2;

                interpolate_to_cell_centre_grid(out u_interp, out v_interp, out w_interp, out p_interp, out div_interp);
            }
            else
            {
                Nx = omega.Nx - 1;
                Ny = omega.Ny - 1;
                Nz = omega.Nz - 1;

                interpolate_to_vertex_grid(out u_interp, out v_interp, out w_interp, out p_interp, out div_interp);
            }

            using (StreamWriter sw = new StreamWriter(fname))
            {
                sw.WriteLine("# vtk DataFile Version 3.0");
                sw.WriteLine("Fast Fluid Dynamics data\n");
                sw.WriteLine("ASCII");
                sw.WriteLine("DATASET RECTILINEAR_GRID");
                sw.WriteLine("FIELD FieldData 1");
                sw.WriteLine("TIME 1 1 double");
                sw.WriteLine("{0}", time);
                sw.WriteLine("DIMENSIONS {0} {1} {2}", omega.Nx - 1, omega.Ny - 1, omega.Nz - 1);
                sw.WriteLine("X_COORDINATES {0} double", omega.Nx - 1);

                for (int i = 0; i < omega.Nx - 1; i++)
                {
                    sw.WriteLine("{0}", hx * i);
                }
                sw.WriteLine("Y_COORDINATES {0} double", omega.Ny - 1);
                for (int i = 0; i < omega.Ny - 1; i++)
                {
                    sw.WriteLine("{0}", hy * i);
                }
                sw.WriteLine("Z_COORDINATES {0} double", omega.Nz - 1);
                for (int i = 0; i < omega.Nz - 1; i++)
                {
                    sw.WriteLine("{0}", hz * i);
                }

                if (cell_centred)
                {
                    sw.WriteLine("CELL_DATA {0}", Nx * Ny * Nz);
                }
                else
                {
                    sw.WriteLine("POINT_DATA {0}", Nx * Ny * Nz);
                }

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
                for (int k = 0; k < Nz; k++)
                {
                    for (int j = 0; j < Ny; j++)
                    {
                        for (int i = 0; i < Nx; i++)
                        {
                            sw.Write("{0} ", p_interp[i, j, k]);
                        }
                    }
                }

                if (cell_centred)
                {
                    sw.WriteLine("SCALARS div_u double {0}", 1);
                    sw.WriteLine("LOOKUP_TABLE default");
                    for (int k = 0; k < Nz; k++)
                    {
                        for (int j = 0; j < Ny; j++)
                        {
                            for (int i = 0; i < Nx; i++)
                            {
                                sw.Write("{0} ", div_interp[i, j, k]);
                            }
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
        * Interpolates pressures and velocities from fs to a cell centred grid
        * 
        * TO DO: implement in parallel 
        ***********************************************************************************/
        private void interpolate_to_cell_centre_grid(out double[, ,] u_interp, out double[, ,] v_interp, 
                    out double[, ,] w_interp, out double[, ,] p_interp, out double[, ,] div_interp)
        {
            int Nx = omega.Nx - 2;
            int Ny = omega.Ny - 2;
            int Nz = omega.Nz - 2;

            u_interp = new double[Nx, Ny, Nz];
            v_interp = new double[Nx, Ny, Nz];
            w_interp = new double[Nx, Ny, Nz];
            p_interp = new double[Nx, Ny, Nz];
            div_interp = new double[Nx, Ny, Nz];

            for (int i = 1; i < Nx + 1; i++)
            {
                for (int j = 1; j < Ny + 1; j++)
                {
                    for (int k = 1; k < Nz + 1; k++)
                    {
                        if (omega.obstacle_cells[i,j,k] == 0)
                        {
                            double[] velocity_interp = de.get_velocity((i - 0.5) * omega.hx, (j - 0.5) * omega.hy, (k - 0.5) * omega.hz); 
                            u_interp[i - 1, j - 1, k - 1] = velocity_interp[0];
                            v_interp[i - 1, j - 1, k - 1] = velocity_interp[1];
                            w_interp[i - 1, j - 1, k - 1] = velocity_interp[2];

                            p_interp[i - 1, j - 1, k - 1] = de.get_pressure((i - 0.5) * omega.hx, (j - 0.5) * omega.hy, (k - 0.5) * omega.hz);

                            div_interp[i - 1, j - 1, k - 1] = (fs.u[i, j, k] - fs.u[i - 1, j, k]) / omega.hx + (fs.v[i, j, k] - fs.v[i, j - 1, k]) / omega.hy +
                                    (fs.w[i, j, k] - fs.w[i, j, k - 1]) / omega.hz;
                         }
                    }
                }
            }
        }

        /************************************************************************************
         * Interpolates pressures and velocities from fs to a vertex centred grid
         * 
         * TO DO: implement in parallel 
         ***********************************************************************************/
        private void interpolate_to_vertex_grid(out double[, ,] u_interp, out double[, ,] v_interp,
                    out double[, ,] w_interp, out double[, ,] p_interp, out double[, ,] div_interp)
        {
            int Nx = omega.Nx - 1;
            int Ny = omega.Ny - 1;
            int Nz = omega.Nz - 1;

            u_interp = new double[Nx, Ny, Nz];
            v_interp = new double[Nx, Ny, Nz];
            w_interp = new double[Nx, Ny, Nz];
            p_interp = new double[Nx, Ny, Nz];
            div_interp = new double[Nx, Ny, Nz];

            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        double[] velocity_interp = de.get_velocity(i * omega.hx, j * omega.hy, k * omega.hz);
                        u_interp[i, j, k] = velocity_interp[0];
                        v_interp[i, j, k] = velocity_interp[1];
                        w_interp[i, j, k] = velocity_interp[2];

                        p_interp[i, j, k] = de.get_pressure(i * omega.hx, j * omega.hy, k* omega.hz);
                    }
                }
            }
        }
    }
}

