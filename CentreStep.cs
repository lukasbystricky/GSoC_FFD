using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolver
{
    class CentreStep : Domain
    {
        public CentreStep(int Nx, int Ny, int Nz)
         {
            this.Nx = Nx;
            this.Ny = Ny;
            this.Nz = Nz;

            length_x = 2;
            length_y = 1;
            length_z = 1;

            hx = length_x / (Nx - 2);
            hy = length_y / (Ny - 2);
            hz = length_z / (Nz - 2);            

            boundary_cells = new int[Nx, Ny, Nz];
            obstacle_cells = new int[Nx, Ny, Nz];
            boundary_normal_x = new int[Nx, Ny, Nz];
            boundary_normal_y = new int[Nx, Ny, Nz];
            boundary_normal_z = new int[Nx, Ny, Nz];
            boundary_u = new double[Nx - 1, Ny, Nz];
            boundary_v = new double[Nx, Ny - 1, Nz];
            boundary_w = new double[Nx, Ny, Nz - 1];

            outflow_cells = new int[Nx, Ny, Nz];

            //C# default values for int or double arrays are 0, so we only need to set nonzero fields

            //create step
            for (int i = (int) Math.Floor(Nx / 3.0); i < (int) Math.Floor(Nx / 2.0); i++)
            {
                for (int j = 1; j < (int)Math.Floor(Ny / 2.0); j++)
                {
                    for (int k = 1; k < Nz - 1; k++)
                    {
                        obstacle_cells[i, j, k] = 1;
                    }
                }
            }

            for (int j = 1; j < (int)Math.Floor(Ny / 2.0) + 1; j++)
            {
                for (int k = 1; k < Nz - 1; k++)
                {
                    boundary_cells[(int)Math.Floor(Nx / 3.0) - 1, j, k] = 1;
                    boundary_cells[(int)Math.Floor(Nx / 2.0), j, k] = 1;

                    boundary_normal_x[(int)Math.Floor(Nx / 3.0) - 1, j, k] = 1;
                    boundary_normal_x[(int)Math.Floor(Nx / 2.0), j, k] = -1;
                }
            }


            for (int i = (int)Math.Floor(Nx / 3.0) - 1; i < (int)Math.Floor(Nx / 2.0) + 2; i ++)
            {
                for (int k = 1; k < Nz - 1; k++)
                {
                    boundary_cells[i, (int)Math.Floor(Ny / 2.0), k] = 1;

                    boundary_normal_y[i, (int)Math.Floor(Ny / 2.0), k] = -1;
                }
            }

            //x = 0 and x = end boundaries
            for (int j = 0; j < Ny; j++)
            {
                for (int k = 0; k < Nz; k++)
                {
                    if (j > 0 && j < Ny - 1 && k > 0 && k < Nz - 1)
                    {
                        boundary_cells[1, j, k] = 1;
                        boundary_cells[Nx - 2, j, k] = 1;
                        boundary_normal_x[1, j, k] = -1;
                        boundary_normal_x[Nx - 2, j, k] = 1;

                        
                        outflow_cells[Nx - 2, j, k] = 1;
                    }
                    
                    boundary_u[0, j, k] = 1;
                    //boundary_u[Nx - 2, j, k] = 1;

                    obstacle_cells[0, j, k] = 1;
                    obstacle_cells[Nx - 1, j, k] = 1;
                }
            }

            //y = 0 and y = end boundaries
            for (int i = 0; i < Nx; i++)
            {
                for (int k = 0; k < Nz; k++)
                {
                    if (i > 0 && i < Nx - 1 && k > 0 && k < Nz - 1)
                    {
                        boundary_cells[i, 1, k] = 1;
                        boundary_cells[i, Ny - 2, k] = 1;
                        
                        boundary_normal_y[i, 1, k] = -1;
                        boundary_normal_y[i, Ny - 2, k] = 1;
                    }

                    obstacle_cells[i, 0, k] = 1;
                    obstacle_cells[i, Ny - 1, k] = 1;                    
                }
            }

            //z = 0 and z = end boundaries
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    if (i > 0 && i < Nx - 1 && j > 0 && j < Ny - 1)
                    {
                        boundary_cells[i, j, 1] = 1;
                        boundary_cells[i, j, Nz - 2] = 1;

                        boundary_normal_z[i, j, 1] = -1;
                        boundary_normal_z[i, j, Nz - 2] = 1;
                    }

                    obstacle_cells[i, j, 0] = 1;
                    obstacle_cells[i, j, Nz - 1] = 1;
                }
            }

             /*******************************************************************
              * Boundary conditions
              *******************************************************************/
            for (int j = 0; j < boundary_u.GetLength(1); j++)
            {
                for (int k = 0; k < boundary_u.GetLength(2); k++)
                {
                    double y, z;
                    if (j == 0)
                    {
                        y = 0;
                    }
                    else if (j == boundary_u.GetLength(1) - 1)
                    {
                        y = length_y;
                    }
                    else
                    {
                        y = (j - 0.5) * hy;
                    }

                    if (k == 0)
                    {
                        z = 0;
                    }
                    else if (k == boundary_u.GetLength(1) - 1)
                    {
                        z = length_z;
                    }
                    else
                    {
                        z = (k - 0.5) * hz;
                    }

                    boundary_u[boundary_u.GetLength(0) - 1, j, k] = z * (z - length_z) * y * (y - length_y);
                    boundary_u[0, j, k] = z * (z - length_z) * y * (y - length_y); 
                }
            }

        }

         /***********************************************************************
          * Copy constructor
          ***********************************************************************/
         public CentreStep(CavityDomain old)
        {
            Nx = old.Nx;
            Ny = old.Ny;
            Nz = old.Nz;

            hx = old.hx;
            hy = old.hy;
            hz = old.hz;

            length_x = old.length_x;
            length_y = old.length_y;
            length_z = old.length_z;

            boundary_cells = old.boundary_cells;
            obstacle_cells = old.obstacle_cells;
            boundary_normal_x = old.boundary_normal_x;
            boundary_normal_y = old.boundary_normal_y;
            boundary_normal_z = old.boundary_normal_z;
            boundary_u = old.boundary_u;
            boundary_v = old.boundary_v;
            boundary_w = old.boundary_w;
        }
    }
}
