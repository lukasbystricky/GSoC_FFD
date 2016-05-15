using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolver
{
    abstract class Domain
    {
        public int N;
        public double h;

        public int[] boundary_nodes;   //flag to indicate if node is on a boundary
        public int[] obstacle;        //flag to indicate if node is part of an obstacle
        public int[] boundary_normal_x; //flag to indicate if boundary at node is normal to x direction
        public int[] boundary_normal_y; //flag to indicate if boundary at node is normal to y direction
        public int[] boundary_normal_z; //flag to indicate if boundary at node is normal to z direction
        public double[] boundary_u;   //x component of velocity at boundary
        public double[] boundary_v;   //y component of velocity at boundary
        public double[] boundary_w;   //z component of velocity at boundary
    }
}
