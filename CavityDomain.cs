using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolver
{
    public class CavityDomain : Domain
    {
        /***************************************************************************
         * Constructor
         **************************************************************************/
        public CavityDomain(int Nx, int Ny, int Nz)
        {
            this.Nx = Nx;
            this.Ny = Ny;
            this.Nz = Nz;

            length_x = 1;
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

            //x = 0 and x = 1 boundaries
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
                    }

                    obstacle_cells[0, j, k] = 1;
                    obstacle_cells[Nx - 1, j, k] = 1;
                }
            }

            //y = 0 and y = 1 boundaries
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

            //z = 0 and z = 1 boundaries
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

            /**************************************************************************************
                * u boundary conditions
             *************************************************************************************/
            for (int i = 0; i < boundary_u.GetLength(0); i++)
            {
                for (int j = 0; j < boundary_u.GetLength(1); j++)
                {
                    boundary_u[i, j, boundary_u.GetLength(2) - 1] = 1;

                    //boundary_u[i, j, 0] = 1;
                }
            }

            /*for (int i = 0; i < boundary_u.GetLength(0); i++)
            {
                for (int k = 0; k < boundary_u.GetLength(2); k++)
                {
                    boundary_u[i, boundary_u.GetLength(1) - 1, k] = 1;

                    boundary_u[i, 0, k] = 1;
                }
            }

            for (int j = 0; j < boundary_u.GetLength(1); j++)
            {
                for (int k = 0; k < boundary_u.GetLength(2); k++)
                {
                    boundary_u[boundary_u.GetLength(0) - 1, j, k] = 1;

                    boundary_u[0, j, k] = 1;
                }
            }*/

            /**************************************************************************************
                * v boundary conditions
             *************************************************************************************/
            /*for (int i = 0; i < boundary_v.GetLength(0); i++)
            {
                for (int j = 0; j < boundary_v.GetLength(1); j++)
                {
                    boundary_v[i, j, boundary_v.GetLength(2) - 1] = 2;

                    boundary_v[i, j, 0] = 2;
                }
            }

            for (int i = 0; i < boundary_v.GetLength(0); i++)
            {
                for (int k = 0; k < boundary_v.GetLength(2); k++)
                {
                    boundary_v[i, boundary_v.GetLength(1) - 1, k] = 2;

                    boundary_v[i, 0, k] = 2;
                }
            }

            for (int j = 0; j < boundary_v.GetLength(1); j++)
            {
                for (int k = 0; k < boundary_v.GetLength(2); k++)
                {
                    boundary_v[boundary_v.GetLength(0) - 1, j, k] = 2;

                    boundary_v[0, j, k] = 2;
                }
            }*/

            /**************************************************************************************
                * w boundary conditions
             *************************************************************************************/
           /*for (int i = 0; i < boundary_w.GetLength(0); i++)
            {
                for (int j = 0; j < boundary_w.GetLength(1); j++)
                {
                    boundary_w[i, j, boundary_w.GetLength(2) - 1] = 3;

                    boundary_w[i, j, 0] = 3;
                }
            }

            for (int i = 0; i < boundary_w.GetLength(0); i++)
            {
                for (int k = 0; k < boundary_w.GetLength(2); k++)
                {
                    boundary_w[i, boundary_w.GetLength(1) - 1, k] = 3;

                    boundary_w[i, 0, k] = 3;
                }
            }

            for (int j = 0; j < boundary_w.GetLength(1); j++)
            {
                for (int k = 0; k < boundary_w.GetLength(2); k++)
                {
                    boundary_w[boundary_w.GetLength(0) - 1, j, k] = 3;

                    boundary_w[0, j, k] = 3;
                }
            }*/
        } 

       

        /***********************************************************************
         * Copy constructor
         ***********************************************************************/
        public CavityDomain(CavityDomain old)
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
