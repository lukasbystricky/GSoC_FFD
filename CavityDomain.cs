using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolver
{
    class CavityDomain : Domain
    {
        int N;
        double h;
        int[] boundary_nodes;   //flag to indicate if node is on a boundary
        int[] obstacle;        //flag to indicate if node is part of an obstacle
        int[] boundary_normal_x; //flag to indicate if boundary at node is normal to x direction
        int[] boundary_normal_y; //flag to indicate if boundary at node is normal to y direction
        int[] boundary_normal_z; //flag to indicate if boundary at node is normal to z direction
        double[] boundary_u;   //x component of velocity at boundary
        double[] boundary_v;   //y component of velocity at boundary
        double[] boundary_w;   //z component of velocity at boundary

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
    }
}
