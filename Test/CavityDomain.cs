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

            set_ghost_flags();

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
