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
        public void export_vtk(String fname)
        {
            using (StreamWriter sw = new StreamWriter(fname))
            {
                sw.WriteLine("# vtk DataFile Version 3.0");
                sw.WriteLine("Fast Fluid Dynamics data\n");
                sw.WriteLine("ASCII");
                sw.WriteLine("DATASET STRUCTURED_GRID");
                sw.WriteLine("DIMENSIONS {0} {1} {2}", N, N, N);
                sw.WriteLine("POINTS {0} double", Math.Pow(N, 3));
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < N; k++)
                        {
                            sw.WriteLine("{0} {1} {2}", h * i, h * j, h * k);
                        }
                    }
                }

                sw.WriteLine("POINT_DATA {0}", Math.Pow(N, 3));
                sw.WriteLine("VECTORS velocity double");
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < N; k++)
                        {
                            sw.WriteLine("{0} {1} {2}", fs.u[FluidSolver.cell_index(i, j, k, N)], 
                                    fs.v[FluidSolver.cell_index(i, j, k, N)], fs.w[FluidSolver.cell_index(i, j, k, N)]);
                        }
                    }
                }

                sw.WriteLine("SCALARS pressure double {0}", 1);
                sw.WriteLine("LOOKUP_TABLE default");
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < N; k++)
                        {
                            sw.WriteLine("{0}", fs.p[cell_index(i, j, k, N)]);
                        }
                    }
                }

                sw.WriteLine("SCALARS nx int {0}", 1);
                sw.WriteLine("LOOKUP_TABLE default");
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < N; k++)
                        {
                            sw.WriteLine("{0}", omega.boundary_normal_x[cell_index(i, j, k, N)]);
                        }
                    }
                }

                sw.WriteLine("SCALARS ny int {0}", 1);
                sw.WriteLine("LOOKUP_TABLE default");
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < N; k++)
                        {
                            sw.WriteLine("{0}", omega.boundary_normal_y[cell_index(i, j, k, N)]);
                        }
                    }
                }

                sw.WriteLine("SCALARS nz int {0}", 1);
                sw.WriteLine("LOOKUP_TABLE default");
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        for (int k = 0; k < N; k++)
                        {
                            sw.WriteLine("{0}", omega.boundary_normal_z[cell_index(i, j, k, N)]);
                        }
                    }
                }
            }
        }
    }
}
