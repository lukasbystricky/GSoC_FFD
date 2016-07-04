using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolver
{
    public class DataExtractor
    {
        private Domain omega;
        private FluidSolver fs;

        public DataExtractor(Domain omega, FluidSolver fs)
        {
            this.omega = omega;
            this.fs = fs;
        }

        public double get_pressure(double x, double y, double z)
        {
            double x_scaled = x / omega.hx;
            double y_scaled = y / omega.hy;
            double z_scaled = z / omega.hz;

            return Utilities.trilinear_interpolation(x_scaled + 0.5,
                                    y_scaled + 0.5, z_scaled + 0.5, fs.p);
        }

        public double[] get_velocity(double x, double y, double z)
        {
            double[] velocity = new double[3];

            double x_scaled = x / omega.hx;
            double y_scaled = y / omega.hy;
            double z_scaled = z / omega.hz;

            velocity[0] = Utilities.trilinear_interpolation(x_scaled, y_scaled + 0.5, z_scaled + 0.5, fs.u);
            velocity[1] = Utilities.trilinear_interpolation(x_scaled + 0.5, y_scaled, z_scaled + 0.5, fs.v);
            velocity[2] = Utilities.trilinear_interpolation(x_scaled + 0.5, y_scaled + 0.5, z_scaled, fs.w);
             
            return velocity;
        }
    }
}
