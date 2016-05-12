using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolver
{
    class FluidSolver
    {
        double[] u;
        double[] v;
        double[] w;
        double[] p;

        void add_force();
        void diffuse();
        void project();
        void advect();
    }
}
