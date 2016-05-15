using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolver
{
    class CavityDomain : Domain
    {
        /***************************************************************************
         * Constructor
         **************************************************************************/
        public CavityDomain(int N)
        {
            this.N = N;
            h = 1.0 / N;

            boundary_nodes = new int[(int) Math.Pow(N,3)];
            obstacle = new int[(int)Math.Pow(N, 3)];
            boundary_normal_x = new int[(int)Math.Pow(N, 3)];
            boundary_normal_y = new int[(int)Math.Pow(N, 3)];
            boundary_normal_z = new int[(int)Math.Pow(N, 3)];
            boundary_u = new double[(int)Math.Pow(N, 3)];
            boundary_v = new double[(int)Math.Pow(N, 3)];
            boundary_w = new double[(int)Math.Pow(N, 3)];

            //C# default values for int or double arrays are 0, so now we only need to set nonzero fields
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    boundary_nodes[FluidSolver.cell_index(0, i, j)] = 1;
                    boundary_nodes[FluidSolver.cell_index(N, i, j)] = 1;

                    boundary_normal_x[FluidSolver.cell_index(0, i, j)] = -1;
                    boundary_normal_x[FluidSolver.cell_index(0, i, j)] = 1;

                    boundary_nodes[FluidSolver.cell_index(i, 0, j)] = 1;
                    boundary_nodes[FluidSolver.cell_index(i, N, j)] = 1;

                    boundary_normal_y[FluidSolver.cell_index(i, 0, j)] = -1;
                    boundary_normal_y[FluidSolver.cell_index(i, N, j)] = 1;

                    boundary_nodes[FluidSolver.cell_index(i, j, 0)] = 1;
                    boundary_nodes[FluidSolver.cell_index(i, j, N)] = 1;

                    boundary_normal_z[FluidSolver.cell_index(i, j, 0)] = -1;
                    boundary_normal_z[FluidSolver.cell_index(i, j, N)] = 1;

                    boundary_u[FluidSolver.cell_index(i, j, N)] = 1;
                }
            }
        }

        /***********************************************************************
         * Copy constructor
         ***********************************************************************/
        public CavityDomain(CavityDomain old)
        {
            N = old.N;
            h = old.h;

            boundary_nodes = old.boundary_nodes;
            obstacle = old.obstacle;
            boundary_normal_x = old.boundary_normal_x;
            boundary_normal_y = old.boundary_normal_y;
            boundary_normal_z = old.boundary_normal_z;
            boundary_u = old.boundary_u;
            boundary_v = old.boundary_v;
            boundary_w = old.boundary_w;
        }
    }
}
