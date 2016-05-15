using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolver
{
    interface Domain
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
    }
}
