using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolver
{
    public class TestDomain : Domain
    {
        private const double a = 1.25;
        private const double d = 2.25;
        private const double nu = 1;

        /***************************************************************************
         * Constructor
         **************************************************************************/
        public TestDomain(int Nx, int Ny, int Nz)
        {
            this.Nx = Nx;
            this.Ny = Ny;
            this.Nz = Nz;

            length_x = 200;
            length_y = 100;
            length_z = 50;

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

            //set up obstacle
            double x_start = 50;
            double x_end = 100;
            double y_start = 40;
            double y_end = 60;
            double z_start = 0;
            double z_end = 30;

            int i_start = (int) Math.Floor(x_start * (Nx - 2) / length_x);
            int i_end = (int)Math.Floor(x_end * (Nx - 2) / length_x);
            int j_start = (int)Math.Floor(y_start * (Ny - 2) / length_y);
            int j_end = (int)Math.Floor(y_end * (Ny - 2) / length_y);
            int k_start = (int)Math.Floor(z_start * (Nz - 2) / length_z);
            int k_end = (int)Math.Floor(z_end * (Nz - 2) / length_z);

            for (int i = i_start; i < i_end; i++)
            {
                for (int j = j_start; j < j_end; j++)
                {
                    for (int k = k_start; k < k_end; k++)
                    {
                        obstacle_cells[i, j, k] = 1;
                    }
                }
            }

            //set outflow boundaries
            //x = 0 will be inflow, z = 0 will be solid ground, all others will be outflow
            
            //x outflow
            for (int j = 1; j < Ny - 1; j++)
            {
                for (int k = 1; k < Nz - 1; k++)
                {
                    outflow_cells[Nx - 2, j, k] = 1;
                }
            }

            //y outflows
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int k = 1; k < Nz - 1; k++)
                {
                    outflow_cells[i, 1, k] = 1;
                    outflow_cells[i, Ny - 2, k] = 1;
                }
            }

            //z outflow
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    outflow_cells[i, j, Nz - 2] = 1;
                }
            }

            set_ghost_flags();
            set_boundary_flags();
            update_boundary_conditions(0);
        }

        public void update_boundary_conditions(double t)
        {
            /**************************************************************************************
            * inflow boundary conditions
            *************************************************************************************/

            for (int j = 0; j < boundary_u.GetLength(1); j++)
            {
                for (int k = 0; k < boundary_u.GetLength(2); k++)
                {
                    boundary_u[0, j, k] = 4;
                }
            }

        }
        /***********************************************************************
         * Copy constructor
         ***********************************************************************/
        public TestDomain(TestDomain old)
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
