using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace FastFluidSolver
{
    /// <summary>
    /// Solves the Navier-Stokes equations using the Fast Fluid Dynamics method
    /// desribed by Stam in the paper "Stable Fluids". 
    /// 
    /// This implementation uses a staggered grid finite difference method to solve the 
    /// spatial equations and backwards Euler in time. Uses a Jacobi iterative method to solve 
    /// the resulting systems. Supports first or second order semi-Lagranian to resolve
    /// advection term. 
    /// 
    /// List of possible future improvements:
    /// 1. Parallelize code, in particular Jacobi solver
    /// 2. Create lists of obstacle and boundary cells to avoid looping over all cells
    ///     when applying boundary conditions
    /// </summary>
    public class FluidSolver
    {
        public struct solver_struct
        {
            public int max_iter; //maximum number of iterations for Gauss-Seidel solver
            public int min_iter; //minimum number of iterations
            public int backtrace_order;
            public double tol; //maximum relative error for Gauss-Seidel solver
            public bool verbose;
        }

        private solver_struct solver_prams;

        public double[, ,] u; // x component of velocity
        public double[, ,] v; // y component of velocity
        public double[, ,] w; // z component of velocity
        public double[, ,] p; // pressure

        private double[, ,] u_old; //scratch arrays for velocities
        private double[, ,] v_old;
        private double[, ,] w_old;
        private double[, ,] p_old;

        private double dt;  //time step
        public int Nx { get; private set; } //number of points in each x coordinate direction
        public int Ny { get; private set; } //number of points in each y coordinate direction
        public int Nz { get; private set; } //number of points in each z coordinate direction

        private double hx;   //spacing in x coordinate direction
        private double hy;   //spacing in x coordinate direction
        private double hz;   //spacing in x coordinate direction
        private double nu;  //fluid viscosity

        Domain omega;

        /// <summary>
        /// Constructor for the fluid solver class.
        /// </summary>
        /// <param name="omega">domain omega</param>
        /// <param name="dt">time step size</param>
        /// <param name="nu">viscosity</param>
        /// <param name="u0">initial x component of velocity</param>
        /// <param name="v0">initial y component of velocity</param>
        /// <param name="w0">initial z component of velocity</param>
        /// <param name="solver_prams">structure containing solver options</param>
        public FluidSolver(Domain omega, double dt, double nu,  double[, ,] u0, double[, ,] v0, 
                double[, ,] w0, solver_struct solver_prams)
        {
            Nx = omega.Nx;
            Ny = omega.Ny;
            Nz = omega.Nz;

            hx = omega.hx;
            hy = omega.hy;
            hz = omega.hz;

            this.dt = dt;
            this.nu = nu;

            this.omega = omega;
            this.solver_prams = solver_prams;

            p = new double[Nx, Ny, Nz];

            u = new double[u0.GetLength(0), u0.GetLength(1), u0.GetLength(2)];
            v = new double[v0.GetLength(0), v0.GetLength(1), v0.GetLength(2)];
            w = new double[w0.GetLength(0), w0.GetLength(1), w0.GetLength(2)];

            u_old = new double[u.GetLength(0), u.GetLength(1), u.GetLength(2)];
            v_old = new double[v.GetLength(0), v.GetLength(1), v.GetLength(2)];
            w_old = new double[w.GetLength(0), w.GetLength(1), w.GetLength(2)];
            p_old = new double[p.GetLength(0), p.GetLength(1), p.GetLength(2)];

            Array.Copy(u0, 0, u, 0, u0.Length);
            Array.Copy(v0, 0, v, 0, v0.Length);
            Array.Copy(w0, 0, w, 0, w0.Length);

            Array.Copy(u, 0, u_old, 0, u.Length);
            Array.Copy(v, 0, v_old, 0, v.Length);
            Array.Copy(w, 0, w_old, 0, w.Length);
        }

        /// <summary>
        /// Copy constructor
        /// </summary>
        public FluidSolver(FluidSolver old)
        {
            Nx = old.Nx;
            Ny = old.Ny;
            Nz = old.Nz;

            hx = old.hx;
            hy = old.hy;
            hz = old.hz;

            dt = old.dt;
            nu = old.nu;

            u = old.u;
            v = old.v;
            w = old.w;
            p = old.p;

            u_old = old.u_old;
            v_old = old.v_old;
            w_old = old.w_old;
            p_old = old.p_old;

            solver_prams = old.solver_prams;
        }


       /// <summary>
       /// Update velocity by adding forcing tem, these can external forces like gravity or 
       /// buoyancy forces for example.
       /// </summary>
       /// <param name="f">forces to add</param>
       /// <param name="x">velocity component array, one of u, v, w</param>
        void add_force(double[, ,] f, ref double[, ,] x)
        {
            int Sx = x.GetLength(0);
            int Sy = x.GetLength(1);
            int Sz = x.GetLength(2);

            for (int i = 0; i < Sx; i++)
            {
                for (int j = 0; j < Sy; j++)
                {
                    for (int k = 0; k < Sz; k++)
                    {
                        x[i, j, k] += dt*f[i, j, k];
                    }
                }
            }
        }

        /// <summary>
        /// Update velocity/concentration by resolving diffusion term. Solves the diffusion 
        /// equation x_t = L(x) using second order finite difference in space and implicit
        /// Euler in time.
        /// </summary>
        /// <param name="x_old">Old state</param>
        /// <param name="x_new">New state</param>
        /// <param name="grid_type">Type of grid used</param>
        /// <remarks>Since we are using a staggered grid, we have to know what kind of 
        /// grid we are using. The following grid numbers are used throughout the program:
        /// 1: cell centred (pressure, temperature, concentration)
        /// 2: centre of faces normal to x direction (x component of velocity)
        /// 3: centre of faces normal to y direction (y component of velocity)
        /// 4: centre of faces normal to z direction (z component of velocity)</remarks>
        void diffuse(double[, ,] x_old, ref double[, ,] x_new, int grid_type)
        {
            double a = 1 + 2 * nu * dt * (Math.Pow(hx, -2) + Math.Pow(hy, -2) + Math.Pow(hz, -2));
            double[] c = new double[6];

            double[, ,] b = new double[x_old.GetLength(0), x_old.GetLength(1), x_old.GetLength(2)];
            Array.Copy(x_old, 0, b, 0, x_old.Length);

            c[0] = -dt * nu * Math.Pow(hz, -2);
            c[1] = -dt * nu * Math.Pow(hy, -2);
            c[2] = -dt * nu * Math.Pow(hx, -2);
            c[3] = c[2];
            c[4] = c[1];
            c[5] = c[0];

            jacobi_solve(a, c, b, x_old, ref x_new, grid_type);
        }

        /// <summary>
        /// Projection step. Solves a Poisson equation for the pressure L(p) = div(u_old)
        /// using finite difference and then updates the velocities u = u_old - grad(p).
        /// </summary>
        void project()
        {
            double[, ,] div = new double[Nx - 1, Ny - 1, Nz - 1];

            // Calculate div(u_old) using finite differences
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        if (omega.obstacle_cells[i, j, k] == 0) 
                        {
                            div[i, j, k] = ((u[i, j, k] - u[i - 1, j, k]) / hx +
                                   (v[i, j, k] - v[i, j - 1, k]) / hy + (w[i, j, k] - w[i, j, k - 1]) / hz) / dt;                     
                        }
                    }
                }
            }

            double a = -2 * (Math.Pow(hx, -2) + Math.Pow(hy, -2) + Math.Pow(hz, -2));
            double[] c = new double[6];

            c[0] = Math.Pow(hz, -2);
            c[1] = Math.Pow(hy, -2);
            c[2] = Math.Pow(hx, -2);
            c[3] = c[2];
            c[4] = c[1];
            c[5] = c[0];

            double[, ,] p0 = new double[Nx, Ny, Nz]; // Initial guess for pressure
            Array.Copy(p, p0, p.Length);

            jacobi_solve(a, c, div, p0, ref p, 1);

            double[] coordinate = new double[3];

            // Update velocity by subtracting grad(p)
            for (int i = 0; i < u.GetLength(0); i++)
            {
                for (int j = 0; j < u.GetLength(1); j++)
                {
                    for (int k = 0; k < u.GetLength(2); k++)
                    {
                        coordinate[0] = i * hx;
                        coordinate[1] = (j - 0.5) * hy;
                        coordinate[2] = (k - 0.5) * hz;

                        if (Utilities.in_domain(coordinate, omega))
                        {
                            u[i, j, k] -= dt * (p[i + 1, j, k] - p[i, j, k]) / hx;
                        }
                    }
                }
            }

            for (int i = 0; i < v.GetLength(0); i++)
            {
                for (int j = 0; j < v.GetLength(1); j++)
                {
                    for (int k = 0; k < v.GetLength(2); k++)
                    {
                        coordinate[0] = (i - 0.5) * hx;
                        coordinate[1] = j * hy;
                        coordinate[2] = (k - 0.5) * hz;

                        if (Utilities.in_domain(coordinate, omega))
                        {
                            v[i, j, k] -= dt * (p[i, j + 1, k] - p[i, j, k]) / hy;
                        }
                    }
                }
            }

            for (int i = 0; i < w.GetLength(0); i++)
            {
                for (int j = 0; j < w.GetLength(1); j++)
                {
                    for (int k = 0; k < w.GetLength(2); k++)
                    {
                        coordinate[0] = (i - 0.5) * hx;
                        coordinate[1] = (j - 0.5) * hy;
                        coordinate[2] = k * hz;

                        if (Utilities.in_domain(coordinate, omega))
                        {
                            w[i, j, k] -= dt * (p[i, j, k + 1] - p[i, j, k]) / hz;
                        }
                    }
                }
            }

            apply_boundary_conditions();
        }

        /// <summary>
        /// Advection step. Resolve advection term by using a semi-Langrangian backtracer.
        /// </summary>
        /// <param name="x">Updated state</param>
        /// <param name="x0">Original state</param>
        /// <param name="velx">x component of velocity</param>
        /// <param name="vely">y component of velocity</param>
        /// <param name="velz">z component of velocity</param>
        /// <param name="grid_type">Grid type, as described in diffusion method</param>
        void advect(ref double[, ,] x, double[, ,] x0, double[, ,] velx, double[, ,] vely, 
                    double[, ,] velz, int grid_type)
        {
            int Sx = x.GetLength(0);
            int Sy = x.GetLength(1);
            int Sz = x.GetLength(2);

            DataExtractor de = new DataExtractor(omega, this);

            // Loop over every node in x
            for (int i = 1; i < Sx - 1; i++)
            {
                for (int j = 1; j < Sy - 1; j++)
                {
                    for (int k = 1; k < Sz - 1; k++)
                    {
                        double xCoord, yCoord, zCoord;
                        
                        double[] velocity0, velocity1;
                        xCoord = yCoord = zCoord = 0;

                        // Get coordinate of node
                        switch (grid_type)
                        {
                            case 1:
                                xCoord = (i - 0.5) * hx;
                                yCoord = (j - 0.5) * hy;
                                zCoord = (k - 0.5) * hz;

                                break;

                            case 2:
                                xCoord = i * hx;
                                yCoord = (j - 0.5) * hy;
                                zCoord = (k - 0.5) * hz;

                                break;

                            case 3:
                                xCoord = (i - 0.5) * hx;
                                yCoord = j * hy;
                                zCoord = (k - 0.5) * hz;

                                break;

                            case 4:
                                xCoord = (i - 0.5) * hx;
                                yCoord = (j - 0.5) * hy;
                                zCoord = k * hz;

                                break;
                        }

                        double[] coordinate = new double[] {xCoord, yCoord, zCoord};

                        if (Utilities.in_domain(coordinate, omega))
                        {
                            // Find velocity at node
                            velocity0 = de.get_velocity(xCoord, yCoord, zCoord);
                            double[] coordBacktraced = new double[3];

                            switch (solver_prams.backtrace_order)
                            {
                                case 1:

                                    // Perform linear backtrace to find origin of fluid element
                                    coordBacktraced[0] = xCoord - dt * (velocity0[0]);
                                    coordBacktraced[1] = yCoord - dt * (velocity0[1]);
                                    coordBacktraced[2] = zCoord - dt * (velocity0[2]);

                                    if (Utilities.in_domain(coordBacktraced, omega))
                                    {
                                        // Set velocity at node to be velocity at backtraced coordinate
                                        x[i, j, k] = Utilities.trilinear_interpolation(i - (dt / hx) * velocity0[0],
                                                    j - (dt / hy) * velocity0[1], k - (dt / hz) * velocity0[2], x0);
                                    }

                                    break;

                                case 2:

                                    // Perform two step second order backtrace to find origin of fluid element
                                    coordBacktraced[0] = xCoord - (dt/2) * (velocity0[0]);
                                    coordBacktraced[1] = yCoord - (dt/2) * (velocity0[1]);
                                    coordBacktraced[2] = zCoord - (dt/2) * (velocity0[2]);

                                    velocity1 = de.get_velocity(coordBacktraced[0], coordBacktraced[1], coordBacktraced[2]);

                                    coordBacktraced[0] -= (dt / 2) * velocity1[0];
                                    coordBacktraced[1] -= (dt / 2) * velocity1[1];
                                    coordBacktraced[2] -= (dt / 2) * velocity1[2];

                                    velocity1 = de.get_velocity(coordBacktraced[0], coordBacktraced[1], coordBacktraced[2]);

                                    if (Utilities.in_domain(coordBacktraced, omega))
                                    {
                                        // Set velocity at node to be velocity at backtraced coordinate
                                        x[i, j, k] = Utilities.trilinear_interpolation(
                                                        i - (dt / (2 * hx)) * (velocity0[0] + velocity1[0]), 
                                                        j - (dt / (2 * hy)) * (velocity0[1] + velocity1[1]),
                                                        k - (dt / (2 * hz)) * (velocity0[2] + velocity1[2]), x0);
                                    }
                                    break;
                            }
                        }
                    }
                }
            }

            apply_boundary_conditions();
        }

        /// <summary>
        /// Applies the boundary conditions. 
        /// </summary>
        /// <remarks>The staggered grid makes this the longest part of the code.</remarks>
        void apply_boundary_conditions()
        {
            // loop over all cells
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    for (int k = 1; k < Nz - 1; k++)
                    {
                        if (omega.obstacle_cells[i,j,k] == 0)
                        {
                            if (omega.boundary_cells[i, j, k] == 1)
                            {
                                /****************************************************************
                                 * 6 faces, +x, -x, +y, -y, +z, -z
                                 * 
                                 * For velocity normal to face, simply prescribe value, for other 
                                 * velocities prescribe average of a point inside the domain and 
                                 * a ghost point outside the domain
                                 ***************************************************************/
                                if (omega.boundary_normal_x[i, j, k] == -1)//-x face
                                {
                                    if (omega.outflow_boundary_x[i, j, k] == 1)
                                    {
                                        p[i - 1, j, k] = p[i, j, k];

                                        u[i - 1, j, k] = u[i, j, k];
                                        v[i - 1, j, k] = v[i, j, k];
                                        w[i - 1, j, k] = w[i, j, k];
                                    }
                                    else
                                    {
                                        p[i - 1, j, k] = p[i, j, k];

                                        u[i - 1, j, k] = omega.boundary_u[i - 1, j, k];
                                        v[i - 1, j, k] = 2 * omega.boundary_v[i - 1, j, k] - v[i, j, k];
                                        w[i - 1, j, k] = 2 * omega.boundary_w[i - 1, j, k] - w[i, j, k];
                                    }
                                }

                                if (omega.boundary_normal_x[i, j, k] == 1)//+x face
                                {
                                    if (omega.outflow_boundary_x[i, j, k] == 1)
                                    {
                                        p[i + 1, j, k] = p[i, j, k];

                                        u[i, j, k] = u[i - 1, j, k];
                                        v[i + 1, j, k] = v[i, j, k];
                                        w[i + 1, j, k] = w[i, j, k];
                                    }
                                    else
                                    {
                                        p[i + 1, j, k] = p[i, j, k];

                                        u[i, j, k] = omega.boundary_u[i, j, k];
                                        v[i + 1, j, k] = 2 * omega.boundary_v[i + 1, j, k] - v[i, j, k];
                                        w[i + 1, j, k] = 2 * omega.boundary_w[i + 1, j, k] - w[i, j, k];
                                    }
                                }

                                if (omega.boundary_normal_y[i, j, k] == -1)//-y face
                                {
                                    if (omega.outflow_boundary_y[i, j, k] == 1)
                                    {
                                        p[i, j - 1, k] = p[i, j, k];

                                        u[i, j - 1, k] = u[i, j, k];
                                        v[i, j - 1, k] = v[i, j, k];
                                        w[i, j - 1, k] = w[i, j, k];
                                    }
                                    else
                                    {
                                        p[i, j - 1, k] = p[i, j, k];

                                        u[i, j - 1, k] = 2 * omega.boundary_u[i, j - 1, k] - u[i, j, k];
                                        v[i, j - 1, k] = omega.boundary_v[i, j - 1, k];
                                        w[i, j - 1, k] = 2 * omega.boundary_w[i, j - 1, k] - w[i, j, k];
                                    }
                                }

                                if (omega.boundary_normal_y[i, j, k] == 1)//+y face
                                {
                                    if (omega.outflow_boundary_y[i, j, k] == 1)
                                    {
                                        p[i, j + 1, k] = p[i, j, k];

                                        u[i, j + 1, k] = u[i, j, k];
                                        v[i, j, k] = v[i, j - 1, k];
                                        w[i, j + 1, k] = w[i, j, k];
                                    }
                                    else
                                    {
                                        p[i, j + 1, k] = p[i, j, k];

                                        u[i, j + 1, k] = 2 * omega.boundary_u[i, j + 1, k] - u[i, j, k];
                                        v[i, j, k] = omega.boundary_v[i, j, k];
                                        w[i, j + 1, k] = 2 * omega.boundary_w[i, j + 1, k] - w[i, j, k];
                                    }
                                }

                                if (omega.boundary_normal_z[i, j, k] == -1)//-z face
                                {
                                    if (omega.outflow_boundary_z[i, j, k] == 1)
                                    {
                                        p[i, j, k - 1] = p[i, j, k];

                                        u[i, j, k - 1] = u[i, j, k];
                                        v[i, j, k - 1] = v[i, j, k];
                                        w[i, j, k - 1] = w[i, j, k];
                                    }
                                    else
                                    {
                                        p[i, j, k - 1] = p[i, j, k];

                                        u[i, j, k - 1] = 2 * omega.boundary_u[i, j, k - 1] - u[i, j, k];
                                        v[i, j, k - 1] = 2 * omega.boundary_v[i, j, k - 1] - v[i, j, k];
                                        w[i, j, k - 1] = omega.boundary_w[i, j, k - 1];
                                    }
                                }

                                if (omega.boundary_normal_z[i, j, k] == 1)//+z face
                                {
                                    if (omega.outflow_boundary_z[i, j, k] == 1)
                                    {
                                        p[i, j, k + 1] = p[i, j, k];

                                        u[i, j, k + 1] = u[i, j, k];
                                        v[i, j, k + 1] = v[i, j, k];
                                        w[i, j, k] = w[i, j, k - 1];
                                    }
                                    else
                                    {
                                        p[i, j, k + 1] = p[i, j, k];

                                        u[i, j, k + 1] = 2 * omega.boundary_u[i, j, k + 1] - u[i, j, k];
                                        v[i, j, k + 1] = 2 * omega.boundary_v[i, j, k + 1] - v[i, j, k];
                                        w[i, j, k] = omega.boundary_w[i, j, k];
                                    }
                                }

                                /********************************************************************
                                 * 12 edges
                                 * 
                                 * For velocities normal to a face, but on an edge where that hasn't
                                 * been assigned yet, prescribe velocity. For velocities tangential to
                                 * the edge, prescribe an average of 4 points around the edge
                                 *******************************************************************/

                                //-x face
                                if (omega.boundary_normal_x[i, j, k] == -1 &&
                                    omega.boundary_normal_y[i, j, k] == -1)
                                {
                                    p[i - 1, j - 1, k] = p[i, j, k];

                                    u[i - 1, j - 1, k] = 2 * omega.boundary_u[i - 1, j - 1, k] -
                                        omega.boundary_u[i - 1, j, k];

                                    v[i - 1, j - 1, k] = 2 * omega.boundary_v[i - 1, j - 1, k] -
                                            omega.boundary_v[i, j - 1, k];

                                    w[i - 1, j - 1, k] = 4 * omega.boundary_w[i - 1, j - 1, k] +
                                        w[i, j, k] - 2 * omega.boundary_w[i, j - 1, k] -
                                        2 * omega.boundary_w[i - 1, j, k];                                    
                                }

                                if (omega.boundary_normal_x[i, j, k] == -1 &&
                                    omega.boundary_normal_y[i, j, k] == 1)
                                {
                                    p[i - 1, j + 1, k] = p[i, j, k];

                                    u[i - 1, j + 1, k] = 2 * omega.boundary_u[i - 1, j + 1, k] -
                                        omega.boundary_u[i - 1, j, k];

                                    v[i - 1, j, k] = 2 * omega.boundary_v[i - 1, j, k] -
                                        omega.boundary_v[i, j, k];

                                    w[i - 1, j + 1, k] = 4 * omega.boundary_w[i - 1, j + 1, k] +
                                        w[i, j, k] - 2 * omega.boundary_w[i - 1, j, k] -
                                        2 * omega.boundary_w[i, j + 1, k];
                                }

                                if (omega.boundary_normal_x[i, j, k] == -1 &&
                                    omega.boundary_normal_z[i, j, k] == -1)
                                {
                                    p[i - 1, j, k - 1] = p[i, j, k];

                                    u[i - 1, j, k - 1] = 2 * omega.boundary_u[i - 1, j, k - 1] -
                                        omega.boundary_u[i - 1, j, k];

                                    w[i - 1, j, k - 1] = 2 * omega.boundary_w[i - 1, j, k - 1] -
                                        omega.boundary_w[i, j, k - 1];

                                    v[i - 1, j, k - 1] = 4 * omega.boundary_v[i - 1, j, k - 1] +
                                        v[i, j, k] - 2 * omega.boundary_v[i - 1, j, k] -
                                        2 * omega.boundary_v[i, j, k - 1];
                                }

                                if (omega.boundary_normal_x[i, j, k] == -1 &&
                                    omega.boundary_normal_z[i, j, k] == 1)
                                {
                                    p[i - 1, j, k + 1] = p[i, j, k];

                                    u[i - 1, j, k + 1] = 2 * omega.boundary_u[i - 1, j, k + 1] -
                                        omega.boundary_u[i - 1, j, k];

                                    v[i - 1, j, k + 1] = 4 * omega.boundary_v[i - 1, j, k + 1] +
                                        v[i, j, k] - 2 * omega.boundary_v[i - 1, j, k] -
                                        2 * omega.boundary_v[i, j, k + 1];

                                    w[i - 1, j, k] = 2 * omega.boundary_w[i - 1, j, k] -
                                        omega.boundary_w[i, j, k];
                                }

                                //+x face
                                if (omega.boundary_normal_x[i, j, k] == 1 &&
                                    omega.boundary_normal_y[i, j, k] == -1)
                                {
                                    p[i + 1, j - 1, k] = p[i, j, k];

                                    u[i, j - 1, k] = 2 * omega.boundary_u[i, j - 1, k] -
                                        omega.boundary_u[i, j, k];

                                    v[i + 1, j - 1, k] = 2 * omega.boundary_v[i + 1, j - 1, k] -
                                        omega.boundary_v[i, j - 1, k];

                                    w[i + 1, j - 1, k] = 4 * omega.boundary_w[i + 1, j - 1, k] +
                                        w[i, j, k] - 2 * omega.boundary_w[i + 1, j, k] -
                                        2 * omega.boundary_w[i, j - 1, k];
                                }

                                if (omega.boundary_normal_x[i, j, k] == 1 &&
                                    omega.boundary_normal_y[i, j, k] == 1)
                                {
                                    p[i + 1, j + 1, k] = p[i, j, k];

                                    u[i, j + 1, k] = 2 * omega.boundary_u[i, j + 1, k] -
                                        omega.boundary_u[i, j, k];

                                    v[i + 1, j, k] = 2 * omega.boundary_v[i + 1, j, k] -
                                        omega.boundary_v[i, j, k];

                                    w[i + 1, j + 1, k] = 4 * omega.boundary_w[i + 1, j + 1, k] +
                                        w[i, j, k] - 2 * omega.boundary_w[i + 1, j, k] -
                                        2 * omega.boundary_w[i, j + 1, k];
                                }

                                if (omega.boundary_normal_x[i, j, k] == 1 &&
                                    omega.boundary_normal_z[i, j, k] == -1)
                                {
                                    p[i + 1, j, k - 1] = p[i, j, k];

                                    u[i, j, k - 1] = 2 * omega.boundary_u[i, j, k - 1] -
                                        omega.boundary_u[i, j, k];

                                    w[i + 1, j, k - 1] = 2 * omega.boundary_w[i + 1, j, k - 1] -
                                        omega.boundary_w[i, j, k - 1];

                                    v[i + 1, j, k - 1] = 4 * omega.boundary_v[i + 1, j, k - 1] +
                                        v[i, j, k] - 2 * omega.boundary_v[i + 1, j, k] -
                                        2 * omega.boundary_v[i, j, k - 1];
                                }

                                if (omega.boundary_normal_x[i, j, k] == 1 &&
                                    omega.boundary_normal_z[i, j, k] == 1)
                                {
                                    p[i + 1, j, k + 1] = p[i, j, k];

                                    u[i, j, k + 1] = 2 * omega.boundary_u[i, j, k + 1] -
                                        omega.boundary_u[i, j, k];

                                    v[i + 1, j, k + 1] = 4 * omega.boundary_v[i + 1, j, k + 1] +
                                        v[i, j, k] - 2 * omega.boundary_v[i + 1, j, k] -
                                        2 * omega.boundary_v[i, j, k + 1];

                                    w[i + 1, j, k] = 2 * omega.boundary_w[i + 1, j, k] -
                                        omega.boundary_w[i, j, k];
                                }

                                //y,z faces
                                if (omega.boundary_normal_y[i, j, k] == -1 &&
                                    omega.boundary_normal_z[i, j, k] == -1)
                                {
                                    p[i, j - 1, k - 1] = p[i, j, k];

                                    u[i, j - 1, k - 1] = 4 * omega.boundary_u[i, j - 1, k - 1] +
                                        u[i, j, k] - 2 * omega.boundary_u[i, j, k - 1] -
                                        2 * omega.boundary_u[i, j - 1, k];

                                    v[i, j - 1, k - 1] = 2 * omega.boundary_v[i, j - 1, k - 1] -
                                        omega.boundary_v[i, j - 1, k];

                                    w[i, j - 1, k - 1] = 2 * omega.boundary_w[i, j - 1, k - 1] -
                                        omega.boundary_w[i, j, k - 1];
                                }

                                if (omega.boundary_normal_y[i, j, k] == -1 &&
                                    omega.boundary_normal_z[i, j, k] == 1)
                                {
                                    p[i, j - 1, k + 1] = p[i, j, k];

                                    u[i, j - 1, k + 1] = 4 * omega.boundary_u[i, j - 1, k + 1] +
                                        u[i, j, k] - 2 * omega.boundary_u[i, j, k + 1] -
                                        2 * omega.boundary_u[i, j - 1, k];

                                    v[i, j - 1, k + 1] = 2 * omega.boundary_v[i, j - 1, k + 1] -
                                        omega.boundary_v[i, j - 1, k];

                                    w[i, j - 1, k] = 2 * omega.boundary_w[i, j - 1, k] -
                                        omega.boundary_w[i, j, k];
                                }

                                if (omega.boundary_normal_y[i, j, k] == 1 &&
                                    omega.boundary_normal_z[i, j, k] == -1)
                                {
                                    p[i, j + 1, k - 1] = p[i, j, k];

                                    u[i, j + 1, k - 1] = 4 * omega.boundary_u[i, j + 1, k - 1] +
                                        u[i, j, k] - 2 * omega.boundary_u[i, j, k - 1] -
                                        2 * omega.boundary_u[i, j + 1, k];

                                    v[i, j, k - 1] = 2 * omega.boundary_v[i, j, k - 1] -
                                        omega.boundary_v[i, j, k];

                                    w[i, j + 1, k - 1] = 2 * omega.boundary_w[i, j + 1, k - 1] -
                                        omega.boundary_w[i, j, k - 1];
                                }

                                if (omega.boundary_normal_y[i, j, k] == 1 &&
                                    omega.boundary_normal_z[i, j, k] == 1)
                                {
                                    p[i, j + 1, k + 1] = p[i, j, k];

                                    u[i, j + 1, k + 1] = 4 * omega.boundary_u[i, j + 1, k + 1] +
                                        u[i, j, k] - 2 * omega.boundary_u[i, j, k + 1] -
                                        2 * omega.boundary_u[i, j + 1, k];

                                    v[i, j, k + 1] = 2 * omega.boundary_v[i, j, k + 1] -
                                        omega.boundary_v[i, j, k];

                                    w[i, j + 1, k] = 2 * omega.boundary_w[i, j + 1, k] -
                                        omega.boundary_w[i, j, k];
                                }

                                /*****************************************************************************
                                 * 8 corners
                                 *****************************************************************************/

                                if (omega.boundary_normal_x[i, j, k] == -1 &&
                                    omega.boundary_normal_y[i, j, k] == -1 &&
                                    omega.boundary_normal_z[i, j, k] == -1)
                                {
                                    u[i - 1, j - 1, k - 1] = 4 * omega.boundary_u[i - 1, j - 1, k - 1] +
                                        u[i - 1, j, k] - 2 * omega.boundary_u[i - 1, j - 1, k] -
                                        2 * omega.boundary_u[i - 1, j, k - 1];

                                    v[i - 1, j - 1, k - 1] = 4 * omega.boundary_v[i - 1, j - 1, k - 1] +
                                        v[i, j - 1, k] - 2 * omega.boundary_v[i - 1, j - 1, k] -
                                        2 * omega.boundary_v[i, j - 1, k - 1];

                                    w[i - 1, j - 1, k - 1] = 4 * omega.boundary_w[i - 1, j - 1, k - 1] +
                                        w[i, j, k - 1] - 2 * omega.boundary_w[i - 1, j, k - 1] -
                                        2 * omega.boundary_w[i, j - 1, k - 1];

                                    p[i - 1, j - 1, k - 1] = p[i, j, k];
                                }

                                if (omega.boundary_normal_x[i, j, k] == -1 &&
                                    omega.boundary_normal_y[i, j, k] == 1 &&
                                    omega.boundary_normal_z[i, j, k] == -1)
                                {
                                    u[i - 1, j + 1, k - 1] = 4 * omega.boundary_u[i - 1, j + 1, k - 1] +
                                        u[i - 1, j, k] - 2 * omega.boundary_u[i - 1, j + 1, k] -
                                        2 * omega.boundary_u[i - 1, j, k - 1];

                                    w[i - 1, j + 1, k - 1] = 4 * omega.boundary_w[i - 1, j + 1, k - 1] +
                                        w[i, j, k - 1] - 2 * omega.boundary_w[i - 1, j, k - 1] -
                                        2 * omega.boundary_w[i, j + 1, k - 1];

                                    p[i - 1, j + 1, k - 1] = p[i, j, k];
                                }

                                if (omega.boundary_normal_x[i, j, k] == -1 &&
                                    omega.boundary_normal_y[i, j, k] == 1 &&
                                    omega.boundary_normal_z[i, j, k] == 1)
                                {
                                    u[i - 1, j + 1, k + 1] = 4 * omega.boundary_u[i - 1, j + 1, k + 1] +
                                        u[i - 1, j, k] - 2 * omega.boundary_u[i - 1, j + 1, k] -
                                        2 * omega.boundary_u[i - 1, j, k + 1];

                                    p[i - 1, j + 1, k + 1] = p[i, j, k];
                                }

                                if (omega.boundary_normal_x[i, j, k] == -1 &&
                                    omega.boundary_normal_y[i, j, k] == -1 &&
                                    omega.boundary_normal_z[i, j, k] == 1)
                                {
                                    u[i - 1, j - 1, k + 1] = 4 * omega.boundary_u[i - 1, j - 1, k + 1] +
                                        u[i - 1, j, k] - 2 * omega.boundary_u[i - 1, j - 1, k] -
                                        2 * omega.boundary_u[i - 1, j, k + 1];

                                    v[i - 1, j - 1, k + 1] = 4 * omega.boundary_v[i - 1, j - 1, k + 1] +
                                        v[i, j - 1, k] - 2 * omega.boundary_v[i - 1, j - 1, k] -
                                        2 * omega.boundary_v[i, j - 1, k + 1];

                                    p[i - 1, j - 1, k + 1] = p[i, j, k];
                                }

                                if (omega.boundary_normal_x[i, j, k] == 1 &&
                                    omega.boundary_normal_y[i, j, k] == -1 &&
                                    omega.boundary_normal_z[i, j, k] == -1)
                                {
                                    v[i + 1, j - 1, k - 1] = 4 * omega.boundary_v[i + 1, j - 1, k - 1] +
                                        v[i, j - 1, k] - 2 * omega.boundary_v[i + 1, j - 1, k] -
                                        2 * omega.boundary_v[i, j - 1, k - 1];

                                    w[i + 1, j - 1, k - 1] = 4 * omega.boundary_w[i + 1, j - 1, k - 1] +
                                        w[i, j, k - 1] - 2 * omega.boundary_w[i + 1, j, k - 1] -
                                        2 * omega.boundary_w[i, j - 1, k - 1];

                                    p[i + 1, j - 1, k - 1] = p[i, j, k];
                                }

                                if (omega.boundary_normal_x[i, j, k] == 1 &&
                                    omega.boundary_normal_y[i, j, k] == 1 &&
                                    omega.boundary_normal_z[i, j, k] == -1)
                                {
                                    w[i + 1, j + 1, k - 1] = 4 * omega.boundary_w[i + 1, j + 1, k - 1] +
                                        w[i, j, k - 1] - 2 * omega.boundary_w[i + 1, j, k - 1] -
                                        2 * omega.boundary_w[i, j + 1, k - 1];

                                    p[i + 1, j + 1, k - 1] = p[i, j, k];
                                }

                                if (omega.boundary_normal_x[i, j, k] == 1 &&
                                    omega.boundary_normal_y[i, j, k] == 1 &&
                                    omega.boundary_normal_z[i, j, k] == 1)
                                {
                                    p[i + 1, j + 1, k + 1] = p[i, j, k];
                                }

                                if (omega.boundary_normal_x[i, j, k] == 1 &&
                                    omega.boundary_normal_y[i, j, k] == -1 &&
                                    omega.boundary_normal_z[i, j, k] == 1)
                                {
                                    v[i + 1, j - 1, k + 1] = 4 * omega.boundary_v[i + 1, j - 1, k + 1] +
                                        v[i, j - 1, k] - 2 * omega.boundary_v[i + 1, j - 1, k] -
                                        2 * omega.boundary_v[i, j - 1, k + 1];

                                    p[i + 1, j - 1, k + 1] = p[i, j, k];
                                }
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Calculate mass inflow and outflow
        /// </summary>
        /// <param name="m_in">Total mass inflow</param>
        /// <param name="m_out">Total mass outflow</param>
        private void calculate_mass_flux(out double m_in, out double m_out)
        {
            m_in = 0;
            m_out = 0;

            // loop over all cells
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    for (int k = 1; k < Nz - 1; k++)
                    {
                        if (omega.obstacle_cells[i, j, k] == 0)
                        {
                            if (omega.boundary_normal_x[i, j, k] == -1)//-x face
                            {
                                if (u[i - 1, j, k] > 0)
                                {
                                    m_in += u[i - 1, j, k] * hy * hz;
                                }
                                else
                                {
                                    m_out -= u[i - 1, j, k] * hy * hz;
                                }
                            }

                            if (omega.boundary_normal_x[i, j, k] == 1)
                            {
                                if (u[i, j, k] > 0)
                                {
                                    m_out += u[i, j, k] * hy * hz;
                                }
                                else
                                {
                                    m_in -= u[i, j, k] * hy * hz;
                                }
                            }

                            if (omega.boundary_normal_y[i, j, k] == -1)
                            {
                                if (v[i, j - 1, k] > 0)
                                {
                                    m_in += v[i, j - 1, k] * hx * hz;
                                }
                                else
                                {
                                    m_out -= v[i, j - 1, k] * hx * hz;
                                }
                            }

                            if (omega.boundary_normal_y[i, j, k] == 1)
                            {
                                if (v[i, j, k] > 0)
                                {
                                    m_out += v[i, j, k] * hx * hz;
                                }
                                else
                                {
                                    m_in -= v[i, j, k] * hx * hz;
                                }
                            }

                            if (omega.boundary_normal_z[i, j, k] == -1)
                            {
                                if (w[i, j - 1, k] > 0)
                                {
                                    m_in += w[i, j, k - 1] * hx * hy;
                                }
                                else
                                {
                                    m_out -= w[i, j, k - 1] * hx * hy;
                                }
                            }

                            if (omega.boundary_normal_z[i, j, k] == 1)
                            {
                                if (w[i, j, k] > 0)
                                {
                                    m_out += w[i, j, k] * hx * hy;
                                }
                                else
                                {
                                    m_in -= w[i, j, k] * hx * hy;
                                }
                            }                            
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Apply mass correction to outflow as described by Zuo et al in 2010 paper
        /// "Improvements in FFD Modeling by Using Different Numerical Schemes"
        /// </summary>
        /// <param name="alpha">Ratio used in correction factor</param>
        /// <remarks>alpha = 1 ensures perfect mass conservation (over domain boundary), 
        /// however paper suggests alpha = 0.7 to avoid instability</remarks>
        private void apply_mass_correction(double alpha)
        {
            double m_in, m_out, correction_factor;
            calculate_mass_flux(out m_in, out m_out);

            correction_factor = 1 + alpha * (m_in / m_out);

            for (int i = 1; i < Nx - 1; i++)
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    for (int k = 1; k < Nz - 1; k++)
                    {
                        if (omega.obstacle_cells[i, j, k] == 0)
                        {
                            if (omega.boundary_normal_x[i, j, k] == -1)
                            {
                                u[i - 1, j, k] *= ((u[i - 1, j, k] < 0) ? correction_factor : 1);  
                            }

                            if (omega.boundary_normal_x[i, j, k] == 1)
                            {
                                u[i, j, k] *= ((u[i, j, k] > 0) ? correction_factor : 1);
                            }

                            if (omega.boundary_normal_y[i, j, k] == -1)
                            {
                                v[i, j - 1, k] *= ((v[i, j - 1, k] < 0) ? correction_factor : 1);
                            }

                            if (omega.boundary_normal_y[i, j, k] == 1)
                            {
                                v[i, j, k] *= ((v[i, j, k] > 0) ? correction_factor : 1);
                            }

                            if (omega.boundary_normal_z[i, j, k] == -1)
                            {
                                w[i, j, k - 1] *= ((w[i, j, k - 1] < 0) ? correction_factor : 1);
                            }

                            if (omega.boundary_normal_x[i, j, k] == 1)
                            {
                                w[i, j, k] *= ((w[i, j, k] > 0) ? correction_factor : 1);
                            }
                        }
                    }
                }
            }
        }


        /// <summary>
        /// Perform a single time step. Add forces, diffuse, project, advect, project.
        /// </summary>
        /// <param name="f_x">x component of forcing term</param>
        /// <param name="f_y">y component of forcing term</param>
        /// <param name="f_z">z component of forcing term</param>
        public void time_step(double[, ,] f_x, double[, ,] f_y, double[, ,] f_z)
        {
            add_force(f_x, ref u);
            add_force(f_y, ref v);
            add_force(f_z, ref w);

            Array.Copy(u, 0, u_old, 0, u.Length);
            Array.Copy(v, 0, v_old, 0, v.Length);
            Array.Copy(w, 0, w_old, 0, w.Length);
            Array.Copy(p, 0, p_old, 0, p.Length);

            diffuse(u_old, ref u, 2);
            diffuse(v_old, ref v, 3);
            diffuse(w_old, ref w, 4);

            project();

            Array.Copy(u, 0, u_old, 0, u.Length);
            Array.Copy(v, 0, v_old, 0, v.Length);
            Array.Copy(w, 0, w_old, 0, w.Length);
            Array.Copy(p, 0, p_old, 0, p.Length);

            advect(ref u, u_old, u_old, v_old, w_old, 2);
            advect(ref v, v_old, u_old, v_old, w_old, 3);
            advect(ref w, w_old, u_old, v_old, w_old, 4);

            apply_mass_correction(0);

            project();
        }

        /// <summary>
        /// Solves the sparse banded system given by the finite difference method applied
        /// to the Poisson or diffusion equation using the iterative Jacobi method.
        /// </summary>
        /// <param name="a">coefficient for diagonal entry</param>
        /// <param name="c">coefficint array other 6 non-zero entries in each row</param>
        /// <param name="b">right hand side</param>
        /// <param name="x0">initial guess</param>
        /// <param name="x1">solution</param>
        /// <param name="grid_type">grid type as described in diffusion method</param>
        /// <remarks>The coefficients for the 6 nonzero entries in each row are given in the 
        /// order x[i,j,k-1], x[i,j-1,k], x[i-1,j,k], x[i+1,j,k], x[i,j+1,k, x[i,j,k+1]</remarks>
        void jacobi_solve(double a, double[] c, double[, ,] b, double[, ,] x0, ref double[, ,] x1, int grid_type)
        {
            int Sx = x0.GetLength(0);
            int Sy = x0.GetLength(1);
            int Sz = x0.GetLength(2);

            int iter = 0;
            double res = 2 * solver_prams.tol;

            double[] coordinate = new double[3];

            while (iter < solver_prams.min_iter || 
                (iter < solver_prams.max_iter && res > solver_prams.tol))
            {
                apply_boundary_conditions();

                for (int k = 1; k < Sz - 1; k++)
                {
                    for (int j = 1; j < Sy - 1; j++)
                    {
                        for (int i = 1; i < Sx - 1; i++)
                        {
                            switch (grid_type)
                            {
                                case 1:
                                    coordinate[0] = (i - 0.5) * hx;
                                    coordinate[1] = (j - 0.5) * hy;
                                    coordinate[2] = (k - 0.5) * hz;
                                    break;

                                case 2:
                                    coordinate[0] = i * hx;
                                    coordinate[1] = (j - 0.5) * hy;
                                    coordinate[2] = (k - 0.5) * hz;
                                    break;

                                case 3:
                                    coordinate[0] = (i - 0.5) * hx;
                                    coordinate[1] = j * hy;
                                    coordinate[2] = (k - 0.5) * hz;
                                    break;

                                case 4:
                                    coordinate[0] = (i - 0.5) * hx;
                                    coordinate[1] = (j - 0.5) * hy;
                                    coordinate[2] = k * hz;
                                    break;
                            }

                            if ( Utilities.in_domain(coordinate, omega)) 
                            {
                                x1[i, j, k] = (b[i, j, k] - (c[0] * x0[i, j, k - 1] +
                                    c[1] * x0[i, j - 1, k] + c[2] * x0[i - 1, j, k] +
                                    c[3] * x0[i + 1, j, k] + c[4] * x0[i, j + 1, k] +
                                    c[5] * x0[i, j, k + 1])) / a;
                                                          
                            }
                        }
                    }
                }

                res = Utilities.compute_L2_difference(x0, x1);
                iter++;

                Array.Copy(x1, 0, x0, 0, x1.Length);
            }     

            if (solver_prams.verbose)
            {
                Console.WriteLine("Jacobi solver completed with residual of {0} in {1} iterations", 
                    res, iter);
            }
        }
    }
}
