using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/*
 * DataExtractor.cs
 * Copyright 2016 Lukas Bystricky <lb13f@my.fsu.edu>
 *
 * This work is licensed under the GNU GPL license version 2 or later.
 */
 
namespace FastFluidSolver
{
    /// <summary>
    /// Extracts velocity and pressure from an FFD simulation
    /// </summary>
    public class DataExtractor
    {
        private Domain omega;
        private FluidSolver fs;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="omega">Domain</param>
        /// <param name="fs">FFD simulation</param>
        public DataExtractor(Domain omega, FluidSolver fs)
        {
            this.omega = omega;
            this.fs = fs;
        }

        /// <summary>
        /// Calculate the pressure at a point (x,y,z)
        /// </summary>
        /// <param name="x">x coordinate</param>
        /// <param name="y">x coordinate</param>
        /// <param name="z">x coordinate</param>
        /// <returns>pressure</returns>
        public double get_pressure(double x, double y, double z)
        {
            double x_scaled = x / omega.hx;
            double y_scaled = y / omega.hy;
            double z_scaled = z / omega.hz;

            return Utilities.trilinear_interpolation(x_scaled + 0.5,
                                    y_scaled + 0.5, z_scaled + 0.5, fs.p);
        }

        /// <summary>
        /// Calculate the velocity (u,v,w) at a point (x,y,z)
        /// </summary>
        /// <param name="x">x coordinate</param>
        /// <param name="y">x coordinate</param>
        /// <param name="z">x coordinate</param>
        /// <returns>velocity (u,v,w)</returns>
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
