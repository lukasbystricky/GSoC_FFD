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

            hx = 1.0 / Nz;
            hy = 1.0 / Nz;
            hz = 1.0 / Nz;

            length_x = 1;
            length_y = 1;
            length_z = 1;

            boundary_cells = new int[Nx, Ny, Nz];
            obstacle_cells = new int[Nx, Ny, Nz];
            boundary_normal_x = new int[Nx, Ny, Nz];
            boundary_normal_y = new int[Nx, Ny, Nz];
            boundary_normal_z = new int[Nx, Ny, Nz];
            boundary_u = new double[Nx, Ny, Nz];
            boundary_v = new double[Nx, Ny, Nz];
            boundary_w = new double[Nx, Ny, Nz];

            //C# default values for int or double arrays are 0, so we only need to set nonzero fields

            //z = 0 and z = 1 boundaries
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    boundary_cells[i, j, 1] = 1;
                    boundary_cells[i, j, Nz - 2] = 1;

                    obstacle_cells[i, j, 0] = 1;
                    obstacle_cells[i, j, Nz - 1] = 1;

                    boundary_normal_z[i, j, 1] = -1;
                    boundary_normal_z[i, j, Nz - 2] = 1;

                    boundary_u[i, j, Nz - 1] = 1;
                }
            }

            //y = 0 and y = 1 boundaries
            for (int i = 0; i < Nx; i++)
            {
                for (int k = 0; k < Nz; k++)
                {
                    boundary_cells[i, 1, k] = 1;
                    boundary_cells[i, Ny - 2, k] = 1;

                    obstacle_cells[i, 0, k] = 1;
                    obstacle_cells[i, Ny - 1, k] = 1;

                    boundary_normal_z[i, 1, k] = -1;
                    boundary_normal_z[i, Ny - 2, k] = 1;
                }
            }

            //x = 0 and x = 1 boundaries
            for (int j = 0; j < Ny; j++)
            {
                for (int k = 0; k < Nz; k++)
                {
                    boundary_cells[1, j, k] = 1;
                    boundary_cells[Nx - 2, j, k] = 1;

                    obstacle_cells[0, j, k] = 1;
                    obstacle_cells[Nx - 1, j, k] = 1;

                    boundary_normal_z[1, j, k] = -1;
                    boundary_normal_z[Nx - 2, j, k] = 1;
                }
            }
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
