using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FastFluidSolver
{
    public abstract class Domain
    {
        public int Nx, Ny, Nz;
        public double hx, hy, hz, length_x, length_y, length_z;

        public int[, ,] boundary_cells { get; protected set; }   //flag to indicate if cell borders a boundary
        public int[, ,] obstacle_cells { get; protected set; }   //flag to indicate if cell is part of an obstacle

        public int[, ,] boundary_normal_x { get; protected set; } //flag to indicate if boundary at cell is normal to x direction
        public int[, ,] boundary_normal_y { get; protected set; } //flag to indicate if boundary at cell is normal to y direction
        public int[, ,] boundary_normal_z { get; protected set; } //flag to indicate if boundary at cell is normal to z direction

        public int[, ,] outflow_cells { get; protected set; }

        public double[, ,] boundary_u { get; protected set; }   //x component of velocity at boundary
        public double[, ,] boundary_v { get; protected set; }   //y component of velocity at boundary
        public double[, ,] boundary_w { get; protected set; }   //z component of velocity at boundary
    }
}
