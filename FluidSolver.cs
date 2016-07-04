using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace FastFluidSolver
{
    /************************************************************************
     * Solves the Navier-Stokes equations using the Fast Fluid Dynamics method
     * desribed by Stam in the paper "Stable Fluids". Uses a staggered grid 
     * finite difference method to solve the spatial equations and backwards
     * Euler in time.
     * 
     * u[i, j, k] located at (i*hx, (j-0.5)*hy, (k-0.5)*hz)
     * i = 0, ..., Nx, j = 0, ... Ny + 1, k = 0, ..., Nz + 1
     * u[0, j, k] and u[Nx-1, j, k] on boundary
     * u[i, 0, k], u[i, j, 0], u[i, Ny+1, k], u[i, j, Nz+1] inside obstacle cell
     * 
     * v[i, j, k] located at ((i-0.5*hx), j*hy, (k-0.5)*hz)
     * i = 0, ..., Nx + 1, j = 0, ... Ny, k = 0, ..., Nz + 1
     * v[i, 0, k] and v[i, Ny-1, k] on boundary
     * v[0, j, k], v[i, j, 0], v[Nx+1, j, k], v[i, j, Nz+1] inside obstacle cell
     * 
     * w[i, j, k] located at((i-0.5)*hx, (j-0.5)*hy, k*hz)
     * i = 0, ..., Nx + 1, j = 0, ... Ny + 1, k = 0, ..., Nz 
     * w[i, j, 0] and w[i, j, Nz-1] on boundary
     * w[0, j, k], w[i, 0, k], w[Nx+1, j, k], w[i, Ny+1, k] inside obstacle cell
     *  
     * p[i,j,k] located at (i*hx, j*hy, k*hz)
     * i = 0, ... Nx + 1, j = 0, ... Ny + 1, k = 0, ... Nz + 1
     * p at i,j,k = 0, Nx+1 are inside obstacle cells
     ************************************************************************/
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
        public double[, ,] v;// y component of velocity
        public double[, ,] w; // z component of velocity
        public double[, ,] p; // pressure

        private double[, ,] u_old; //scratch arrays for velocities
        private double[, ,] v_old;
        private double[, ,] w_old;

        private double dt;  //time step
        public int Nx { get; private set; } //number of points in each x coordinate direction
        public int Ny { get; private set; } //number of points in each y coordinate direction
        public int Nz { get; private set; } //number of points in each z coordinate direction

        private double hx;   //spacing in x coordinate direction
        private double hy;   //spacing in x coordinate direction
        private double hz;   //spacing in x coordinate direction
        private double nu;  //fluid viscosity

        Domain omega;

        void initialize() { }
        void add_force() { }

        /****************************************************************************
         * Constructor
         ****************************************************************************/
        public FluidSolver(Domain omega, double dt, double nu, double[, ,] u0, double[, ,] v0, double[, ,] w0, solver_struct solver_prams)
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

            Array.Copy(u0, 0, u, 0, u0.Length);
            Array.Copy(v0, 0, v, 0, v0.Length);
            Array.Copy(w0, 0, w, 0, w0.Length);

            Array.Copy(u, 0, u_old, 0, u.Length);
            Array.Copy(v, 0, v_old, 0, v.Length);
            Array.Copy(w, 0, w_old, 0, w.Length);
        }

        /****************************************************************************
         * Copy constructor
         ****************************************************************************/
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

            solver_prams = old.solver_prams;
        }

        /****************************************************************************
         * Diffusion step. Diffuse solve diffusion equation x_t = L(x) using second 
         * order finite difference in space and backwards Euler in time
         ****************************************************************************/
        void diffuse(double[, ,] x_old, ref double[, ,] x_new)
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

            gs_solve(a, c, b, x_old, ref x_new);
        }

        /*****************************************************************************
         * Projection step. Solves Poisson equation L(p) = div(u_old) using finte
         * difference and updates the velocities, u = u_old - grad(p)
         ****************************************************************************/
        void project()
        {
            double[, ,] div = new double[Nx - 1, Ny - 1, Nz - 1];

            // calculate div(w) using finite differences
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        if (omega.obstacle_cells[i, j, k] == 0) //node not inside an obstacle
                        {
                            div[i, j, k] = (u[i, j, k] - u[i - 1, j, k]) / hx +
                                   (v[i, j, k] - v[i, j - 1, k]) / hy + (w[i, j, k] - w[i, j, k - 1]) / hz;                      
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

            double[, ,] p0 = new double[Nx, Ny, Nz]; //initial guess for pressure
            Array.Copy(p, p0, p.Length);

            gs_solve(a, c, div, p0, ref p);

            //update velocity by subtracting calculate grad(p)
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        if (omega.obstacle_cells[i, j, k] == 0) //node not inside an obstacle
                        {
                            u[i, j, k] -= dt*(p[i + 1, j, k] - p[i, j, k]) / hx;
                            v[i, j, k] -= dt*(p[i, j + 1, k] - p[i, j, k]) / hy;
                            w[i, j, k] -= dt*(p[i, j, k + 1] - p[i, j, k]) / hz;
                        }
                    }
                }
            }

            apply_boundary_conditions();
        }

        /***************************************************************************
         * Advection step. Uses a first order backtrace to update x.
         * 
         * @inputs
         * double[, ,] x - reference to quantity to advect, can be a velocity, a
         * concentration or temperature
         * double[, ,] x0 - initial state of x before advection
         * double[, ,] velx - x velocity
         * double[, ,] vely - y velocity 
         * double[, ,] velz - z velocity 
         * int grid_type - specifies type of grid: 1 for cell centred, 2 for u velocity
         *         3 for v velocity, 4 for w velocity
         *         
         * TO DO: implement in parallel
         **************************************************************************/
        void advect(ref double[, ,] x, double[, ,] x0, double[, ,] velx, double[, ,] vely,  double[, ,] velz, int grid_type)
        {
            int Sx = x.GetLength(0);
            int Sy = x.GetLength(1);
            int Sz = x.GetLength(2);

            DataExtractor de = new DataExtractor(omega, this);

            for (int i = 1; i < Sx - 1; i++)
            {
                for (int j = 1; j < Sy - 1; j++)
                {
                    for (int k = 1; k < Sz - 1; k++)
                    {
                        if (omega.obstacle_cells[i, j, k] == 0)
                        {

                            double xCoord, yCoord, zCoord;
                            double[] velocity0, velocity1;

                            xCoord = yCoord = zCoord = 0;

                            switch (grid_type)
                            {
                                case 1:
                                    xCoord = (i - 0.5) * hx;
                                    yCoord = (j - 0.5) * hy;
                                    zCoord = (k - 0.5) * hz;

                                    break;

                                case 2:
                                    xCoord = (i - 1) * hx;
                                    yCoord = (j - 0.5) * hy;
                                    zCoord = (k - 0.5) * hz;

                                    break;

                                case 3:
                                    xCoord = (i - 0.5) * hx;
                                    yCoord = (j - 1) * hy;
                                    zCoord = (k - 0.5) * hz;

                                    break;

                                case 4:
                                    xCoord = (i - 0.5) * hx;
                                    yCoord = (j - 0.5) * hy;
                                    zCoord = (k - 1) * hz;

                                    break;
                            }

                            velocity0 = de.get_velocity(xCoord, yCoord, zCoord);
                            double[] coordNew = new double[3];

                            switch (solver_prams.backtrace_order)
                            {
                                case 1:

                                    coordNew[0] = xCoord - dt * (velocity0[0]);
                                    coordNew[1] = yCoord - dt * (velocity0[1]);
                                    coordNew[2] = zCoord - dt * (velocity0[2]);

                                    if (Utilities.in_domain(coordNew, omega))
                                    {
                                        x[i, j, k] = Utilities.trilinear_interpolation(i - (dt / hx) * velocity0[0],
                                                    j - (dt / hy) * velocity0[1], k - (dt / hz) * velocity0[2], x0);
                                    }

                                    break;

                                case 2:

                                    velocity1 = de.get_velocity(xCoord - dt * velocity0[0], yCoord - dt * velocity0[1], zCoord - dt * velocity0[2]);

                                    coordNew[0] = xCoord - (dt / 2) * (velocity0[0] - velocity1[0]);
                                    coordNew[1] = yCoord - (dt / 2) * (velocity0[1] - velocity1[1]);
                                    coordNew[2] = zCoord - (dt / 2) * (velocity0[2] - velocity1[2]);

                                    if (Utilities.in_domain(coordNew, omega))
                                    {
                                        x[i, j, k] = Utilities.trilinear_interpolation(i - (dt / (2 * hx)) * (velocity0[0] + velocity1[0]),
                                                    j - (dt / (2 * hy)) * (velocity0[1] + velocity1[1]), k - (dt / (2 * hz)) * (velocity0[2] + velocity1[2]), x0);
                                    }
                                    break;
                            }
                        }                  
                    }
                }
            }

            apply_boundary_conditions();
        }

        /*********************************************************************************
         * Applies the Dirichlet boundary conditions from the domain omega to the velocities
         * and homogeneuous Neumann boundary conditions to the pressure
         ********************************************************************************/
        void apply_boundary_conditions()
        {
            //loop over all cells
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    for (int k = 0; k < Nz; k++)
                    {
                        if (omega.obstacle_cells[i,j,k] == 0)
                        {
                            if (omega.boundary_cells[i, j, k] == 1)
                            {
                                if (omega.outflow_cells[i, j, k] == 0)
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
                                        p[i - 1, j, k] = p[i, j, k];

                                        u[i - 1, j, k] = omega.boundary_u[i - 1, j, k];
                                        v[i - 1, j, k] = 2 * omega.boundary_v[i - 1, j, k] - v[i, j, k];
                                        w[i - 1, j, k] = 2 * omega.boundary_w[i - 1, j, k] - w[i, j, k];
                                    }

                                    if (omega.boundary_normal_x[i, j, k] == 1)//+x face
                                    {
                                        p[i + 1, j, k] = p[i, j, k];

                                        u[i, j, k] = omega.boundary_u[i, j, k];
                                        v[i + 1, j, k] = 2 * omega.boundary_v[i + 1, j, k] - v[i, j, k];
                                        w[i + 1, j, k] = 2 * omega.boundary_w[i + 1, j, k] - w[i, j, k];
                                    }

                                    if (omega.boundary_normal_y[i, j, k] == -1)//-y face
                                    {
                                        p[i, j - 1, k] = p[i, j, k];

                                        u[i, j - 1, k] = 2 * omega.boundary_u[i, j - 1, k] - u[i, j, k];
                                        v[i, j - 1, k] = omega.boundary_v[i, j - 1, k];
                                        w[i, j - 1, k] = 2 * omega.boundary_w[i, j - 1, k] - w[i, j, k];
                                    }

                                    if (omega.boundary_normal_y[i, j, k] == 1)//+y face
                                    {
                                        p[i, j + 1, k] = p[i, j, k];

                                        u[i, j + 1, k] = 2 * omega.boundary_u[i, j + 1, k] - u[i, j, k];
                                        v[i, j, k] = omega.boundary_v[i, j, k];
                                        w[i, j + 1, k] = 2 * omega.boundary_w[i, j + 1, k] - w[i, j, k];
                                    }

                                    if (omega.boundary_normal_z[i, j, k] == -1)//-z face
                                    {
                                        p[i, j, k - 1] = p[i, j, k];

                                        u[i, j, k - 1] = 2 * omega.boundary_u[i, j, k - 1] - u[i, j, k];
                                        v[i, j, k - 1] = 2 * omega.boundary_v[i, j, k - 1] - v[i, j, k];
                                        w[i, j, k - 1] = omega.boundary_w[i, j, k - 1];
                                    }

                                    if (omega.boundary_normal_z[i, j, k] == 1)//+z face
                                    {
                                        p[i, j, k + 1] = p[i, j, k];

                                        u[i, j, k + 1] = 2 * omega.boundary_u[i, j, k + 1] - u[i, j, k];
                                        v[i, j, k + 1] = 2 * omega.boundary_v[i, j, k + 1] - v[i, j, k];
                                        w[i, j, k] = omega.boundary_w[i, j, k];
                                    }

                                    /********************************************************************
                                     * 12 edges
                                     * 
                                     * For velocities normal to a face, but on an edge where that hasn't
                                     * been assigned yet, prescribe velocity. For velocities tangential to
                                     * the edge, prescribe an average of 4 points around the edge
                                     *******************************************************************/

                                    //-x face
                                    if (omega.boundary_normal_x[i, j, k] == -1 && omega.boundary_normal_y[i, j, k] == -1)
                                    {
                                        p[i - 1, j - 1, k] = p[i, j, k];

                                        u[i - 1, j - 1, k] = 2 * omega.boundary_u[i - 1, j - 1, k] - omega.boundary_u[i - 1, j, k];
                                        v[i - 1, j - 1, k] = 2 * omega.boundary_v[i - 1, j - 1, k] - omega.boundary_v[i, j - 1, k];
                                        w[i - 1, j - 1, k] = 4 * omega.boundary_w[i - 1, j - 1, k] + w[i, j, k] - 2 * omega.boundary_w[i, j - 1, k] -
                                                                    2 * omega.boundary_w[i - 1, j, k];
                                    }

                                    if (omega.boundary_normal_x[i, j, k] == -1 && omega.boundary_normal_y[i, j, k] == 1)
                                    {
                                        p[i - 1, j + 1, k] = p[i, j, k];

                                        u[i - 1, j + 1, k] = 2 * omega.boundary_u[i - 1, j + 1, k] - omega.boundary_u[i - 1, j, k];
                                        v[i - 1, j, k] = 2 * omega.boundary_v[i - 1, j, k] - omega.boundary_v[i, j, k];
                                        w[i - 1, j + 1, k] = 4 * omega.boundary_w[i - 1, j + 1, k] + w[i, j, k] - 2 * omega.boundary_w[i - 1, j, k] -
                                                                    2 * omega.boundary_w[i, j + 1, k];
                                    }

                                    if (omega.boundary_normal_x[i, j, k] == -1 && omega.boundary_normal_z[i, j, k] == -1)
                                    {
                                        p[i - 1, j, k - 1] = p[i, j, k];

                                        u[i - 1, j, k - 1] = 2 * omega.boundary_u[i - 1, j, k - 1] - omega.boundary_u[i - 1, j, k];
                                        w[i - 1, j, k - 1] = 2 * omega.boundary_w[i - 1, j, k - 1] - omega.boundary_w[i, j, k - 1];
                                        v[i - 1, j, k - 1] = 4 * omega.boundary_v[i - 1, j, k - 1] + v[i, j, k] - 2 * omega.boundary_v[i - 1, j, k] -
                                                                        2 * omega.boundary_v[i, j, k - 1];
                                    }

                                    if (omega.boundary_normal_x[i, j, k] == -1 && omega.boundary_normal_z[i, j, k] == 1)
                                    {
                                        p[i - 1, j, k + 1] = p[i, j, k];

                                        u[i - 1, j, k + 1] = 2 * omega.boundary_u[i - 1, j, k + 1] - omega.boundary_u[i - 1, j, k];
                                        v[i - 1, j, k + 1] = 4 * omega.boundary_v[i - 1, j, k + 1] + v[i, j, k] - 2 * omega.boundary_v[i - 1, j, k] -
                                                                        2 * omega.boundary_v[i, j, k + 1];
                                        w[i - 1, j, k] = 2 * omega.boundary_w[i - 1, j, k] - 2 * omega.boundary_w[i, j, k];
                                    }

                                    //+x face
                                    if (omega.boundary_normal_x[i, j, k] == 1 && omega.boundary_normal_y[i, j, k] == -1)
                                    {
                                        p[i + 1, j - 1, k] = p[i, j, k];

                                        u[i, j - 1, k] = 2 * omega.boundary_u[i, j - 1, k] - omega.boundary_u[i, j, k];
                                        v[i + 1, j - 1, k] = 2 * omega.boundary_v[i + 1, j - 1, k] - omega.boundary_v[i, j - 1, k];
                                        w[i + 1, j - 1, k] = 4 * omega.boundary_w[i + 1, j - 1, k] + w[i, j, k] - 2 * omega.boundary_w[i + 1, j, k] -
                                                                        2 * omega.boundary_w[i, j - 1, k];
                                    }

                                    if (omega.boundary_normal_x[i, j, k] == 1 && omega.boundary_normal_y[i, j, k] == 1)
                                    {
                                        p[i + 1, j + 1, k] = p[i, j, k];

                                        u[i, j + 1, k] = 2 * omega.boundary_u[i, j + 1, k] - omega.boundary_u[i, j, k];
                                        v[i + 1, j, k] = 2 * omega.boundary_v[i + 1, j, k] - omega.boundary_v[i, j, k];
                                        w[i + 1, j + 1, k] = 4 * omega.boundary_w[i + 1, j + 1, k] + w[i, j, k] - 2 * omega.boundary_w[i + 1, j, k] -
                                                                        2 * omega.boundary_w[i, j + 1, k];
                                    }

                                    if (omega.boundary_normal_x[i, j, k] == 1 && omega.boundary_normal_z[i, j, k] == -1)
                                    {
                                        p[i + 1, j, k - 1] = p[i, j, k];

                                        u[i, j, k - 1] = 2 * omega.boundary_u[i, j, k - 1] - omega.boundary_u[i, j, k];
                                        w[i + 1, j, k - 1] = 2 * omega.boundary_w[i + 1, j, k - 1] - omega.boundary_w[i, j, k - 1];
                                        v[i + 1, j, k - 1] = 4 * omega.boundary_v[i + 1, j, k - 1] + v[i, j, k] - 2 * omega.boundary_v[i + 1, j, k] -
                                                                        2 * omega.boundary_v[i, j, k - 1];
                                    }

                                    if (omega.boundary_normal_x[i, j, k] == 1 && omega.boundary_normal_z[i, j, k] == 1)
                                    {
                                        p[i + 1, j, k + 1] = p[i, j, k];

                                        u[i, j, k + 1] = 2 * omega.boundary_u[i, j, k + 1] - omega.boundary_u[i, j, k];
                                        v[i + 1, j, k + 1] = 4 * omega.boundary_v[i + 1, j, k + 1] + v[i, j, k] - 2 * omega.boundary_v[i + 1, j, k] -
                                                                        2 * omega.boundary_v[i, j, k + 1];
                                        w[i + 1, j, k] = 2 * omega.boundary_w[i + 1, j, k] - omega.boundary_w[i, j, k];
                                    }

                                    //y,z faces
                                    if (omega.boundary_normal_y[i, j, k] == -1 && omega.boundary_normal_z[i, j, k] == -1)
                                    {
                                        p[i, j - 1, k - 1] = p[i, j, k];

                                        u[i, j - 1, k - 1] = 4 * omega.boundary_u[i, j - 1, k - 1] + u[i, j, k] - 2 * omega.boundary_u[i, j, k - 1] -
                                                                    2 * omega.boundary_u[i, j - 1, k];
                                        v[i, j - 1, k - 1] = 2 * omega.boundary_v[i, j - 1, k - 1] - omega.boundary_v[i, j - 1, k];
                                        w[i, j - 1, k - 1] = 2 * omega.boundary_w[i, j - 1, k - 1] - omega.boundary_w[i, j, k - 1];
                                    }

                                    if (omega.boundary_normal_y[i, j, k] == -1 && omega.boundary_normal_z[i, j, k] == 1)
                                    {
                                        p[i, j - 1, k + 1] = p[i, j, k];

                                        u[i, j - 1, k + 1] = 4 * omega.boundary_u[i, j - 1, k + 1] + u[i, j, k] - 2 * omega.boundary_u[i, j, k + 1] -
                                                                    2 * omega.boundary_u[i, j - 1, k];
                                        v[i, j - 1, k + 1] = 2 * omega.boundary_v[i, j - 1, k + 1] - omega.boundary_v[i, j - 1, k];
                                        w[i, j - 1, k] = 2 * omega.boundary_w[i, j - 1, k] - omega.boundary_w[i, j, k];
                                    }

                                    if (omega.boundary_normal_y[i, j, k] == 1 && omega.boundary_normal_z[i, j, k] == -1)
                                    {
                                        p[i, j + 1, k - 1] = p[i, j, k];

                                        u[i, j + 1, k - 1] = 4 * omega.boundary_u[i, j + 1, k - 1] + u[i, j, k] - 2 * omega.boundary_u[i, j, k - 1] -
                                                                    2 * omega.boundary_u[i, j + 1, k];
                                        v[i, j, k - 1] = 2 * omega.boundary_v[i, j, k - 1] - omega.boundary_v[i, j, k];
                                        w[i, j + 1, k - 1] = 2 * omega.boundary_w[i, j + 1, k - 1] - omega.boundary_w[i, j, k - 1];
                                    }

                                    if (omega.boundary_normal_y[i, j, k] == 1 && omega.boundary_normal_z[i, j, k] == 1)
                                    {
                                        p[i, j + 1, k + 1] = p[i, j, k];

                                        u[i, j + 1, k + 1] = 4 * omega.boundary_u[i, j + 1, k + 1] + u[i, j, k] - 2 * omega.boundary_u[i, j, k + 1] -
                                                                    2 * omega.boundary_u[i, j + 1, k];
                                        v[i, j, k + 1] = 2 * omega.boundary_v[i, j, k + 1] - omega.boundary_v[i, j, k];
                                        w[i, j + 1, k] = 2 * omega.boundary_w[i, j + 1, k] - omega.boundary_w[i, j, k];
                                    }

                                    /*****************************************************************************
                                     * 8 corners
                                     *****************************************************************************/

                                    if (omega.boundary_normal_x[i, j, k] == -1 && omega.boundary_normal_y[i, j, k] == -1 && omega.boundary_normal_z[i, j, k] == -1)
                                    {
                                        u[i - 1, j - 1, k - 1] = 4 * omega.boundary_u[i - 1, j - 1, k - 1] + u[i - 1, j, k] - 2 * omega.boundary_u[i - 1, j - 1, k] -
                                                                                 2 * omega.boundary_u[i - 1, j, k - 1];

                                        v[i - 1, j - 1, k - 1] = 4 * omega.boundary_v[i - 1, j - 1, k - 1] + v[i, j - 1, k] - 2 * omega.boundary_v[i - 1, j - 1, k] -
                                                                                 2 * omega.boundary_v[i, j - 1, k - 1];

                                        w[i - 1, j - 1, k - 1] = 4 * omega.boundary_w[i - 1, j - 1, k - 1] + w[i, j, k - 1] - 2 * omega.boundary_w[i - 1, j, k - 1] -
                                                                                 2 * omega.boundary_w[i, j - 1, k - 1];

                                        p[i - 1, j - 1, k - 1] = p[i, j, k];
                                    }

                                    if (omega.boundary_normal_x[i, j, k] == -1 && omega.boundary_normal_y[i, j, k] == 1 && omega.boundary_normal_z[i, j, k] == -1)
                                    {
                                        u[i - 1, j + 1, k - 1] = 4 * omega.boundary_u[i - 1, j + 1, k - 1] + u[i - 1, j, k] - 2 * omega.boundary_u[i - 1, j + 1, k] -
                                                                                 2 * omega.boundary_u[i - 1, j, k - 1];

                                        w[i - 1, j + 1, k - 1] = 4 * omega.boundary_w[i - 1, j + 1, k - 1] + w[i, j, k - 1] - 2 * omega.boundary_w[i - 1, j, k - 1] -
                                                                                 2 * omega.boundary_w[i, j + 1, k - 1];

                                        p[i - 1, j + 1, k - 1] = p[i, j, k];
                                    }

                                    if (omega.boundary_normal_x[i, j, k] == -1 && omega.boundary_normal_y[i, j, k] == 1 && omega.boundary_normal_z[i, j, k] == 1)
                                    {
                                        u[i - 1, j + 1, k + 1] = 4 * omega.boundary_u[i - 1, j + 1, k + 1] + u[i - 1, j, k] - 2 * omega.boundary_u[i - 1, j + 1, k] -
                                                                                2 * omega.boundary_u[i - 1, j, k + 1];
                                        p[i - 1, j + 1, k + 1] = p[i, j, k];
                                    }

                                    if (omega.boundary_normal_x[i, j, k] == -1 && omega.boundary_normal_y[i, j, k] == -1 && omega.boundary_normal_z[i, j, k] == 1)
                                    {
                                        u[i - 1, j - 1, k + 1] = 4 * omega.boundary_u[i - 1, j - 1, k + 1] + u[i - 1, j, k] - 2 * omega.boundary_u[i - 1, j - 1, k] -
                                                                                2 * omega.boundary_u[i - 1, j, k + 1];

                                        v[i - 1, j - 1, k + 1] = 4 * omega.boundary_v[i - 1, j - 1, k + 1] + v[i, j - 1, k] - 2 * omega.boundary_v[i - 1, j - 1, k] -
                                                                                2 * omega.boundary_v[i, j - 1, k + 1];

                                        p[i - 1, j - 1, k + 1] = p[i, j, k];
                                    }

                                    if (omega.boundary_normal_x[i, j, k] == 1 && omega.boundary_normal_y[i, j, k] == -1 && omega.boundary_normal_z[i, j, k] == -1)
                                    {
                                        v[i + 1, j - 1, k - 1] = 4 * omega.boundary_v[i + 1, j - 1, k - 1] + v[i, j - 1, k] - 2 * omega.boundary_v[i + 1, j - 1, k] -
                                                                                 2 * omega.boundary_v[i, j - 1, k - 1];

                                        w[i + 1, j - 1, k - 1] = 4 * omega.boundary_w[i + 1, j - 1, k - 1] + w[i, j, k - 1] - 2 * omega.boundary_w[i + 1, j, k - 1] -
                                                                                 omega.boundary_w[i, j - 1, k - 1];

                                        p[i + 1, j - 1, k - 1] = p[i, j, k];
                                    }

                                    if (omega.boundary_normal_x[i, j, k] == 1 && omega.boundary_normal_y[i, j, k] == 1 && omega.boundary_normal_z[i, j, k] == -1)
                                    {
                                        w[i + 1, j + 1, k - 1] = 4 * omega.boundary_w[i + 1, j + 1, k - 1] + w[i, j, k - 1] - 2 * omega.boundary_w[i + 1, j, k - 1] -
                                                                                 2 * omega.boundary_w[i, j + 1, k - 1];

                                        p[i + 1, j + 1, k - 1] = p[i, j, k];
                                    }

                                    if (omega.boundary_normal_x[i, j, k] == 1 && omega.boundary_normal_y[i, j, k] == 1 && omega.boundary_normal_z[i, j, k] == 1)
                                    {
                                        p[i + 1, j + 1, k + 1] = p[i, j, k];
                                    }

                                    if (omega.boundary_normal_x[i, j, k] == 1 && omega.boundary_normal_y[i, j, k] == -1 && omega.boundary_normal_z[i, j, k] == 1)
                                    {
                                        v[i + 1, j - 1, k + 1] = 4 * omega.boundary_v[i + 1, j - 1, k + 1] + v[i, j - 1, k] - 2 * omega.boundary_v[i + 1, j - 1, k] -
                                                                                 2 * omega.boundary_v[i, j - 1, k + 1];

                                        p[i + 1, j - 1, k + 1] = p[i, j, k];
                                    }
                                }
                                else
                                {
                                    if (omega.boundary_normal_x[i, j, k] == -1)//-x face
                                    {
                                        p[i - 1, j, k] = p[i, j, k];

                                        u[i - 1, j, k] = u[i, j, k];
                                        v[i - 1, j, k] = v[i, j, k];
                                        w[i - 1, j, k] = w[i, j, k];
                                    }

                                    if (omega.boundary_normal_x[i, j, k] == 1)//+x face
                                    {
                                        p[i + 1, j, k] = p[i, j, k];

                                        u[i, j, k] = u[i - 1, j, k];
                                        v[i + 1, j, k] = v[i, j, k];
                                        w[i + 1, j, k] = w[i, j, k];
                                    }

                                    if (omega.boundary_normal_y[i, j, k] == -1)//-y face
                                    {
                                        p[i, j - 1, k] = p[i, j, k];

                                        u[i, j - 1, k] = u[i, j, k];
                                        v[i, j - 1, k] = v[i, j, k];
                                        w[i, j - 1, k] = w[i, j, k];
                                    }

                                    if (omega.boundary_normal_y[i, j, k] == 1)//+y face
                                    {
                                        p[i, j + 1, k] = p[i, j, k];

                                        u[i, j + 1, k] = u[i, j, k];
                                        v[i, j, k] = v[i, j - 1, k];
                                        w[i, j + 1, k] = w[i, j, k];
                                    }

                                    if (omega.boundary_normal_z[i, j, k] == -1)//-z face
                                    {
                                        p[i, j, k - 1] = p[i, j, k];

                                        u[i, j, k - 1] = u[i, j, k];
                                        v[i, j, k - 1] = v[i, j, k];
                                        w[i, j, k - 1] = w[i, j, k];
                                    }

                                    if (omega.boundary_normal_z[i, j, k] == 1)//+z face
                                    {
                                        p[i, j, k + 1] = p[i, j, k];

                                        u[i, j, k + 1] = u[i, j, k];
                                        v[i, j, k + 1] = v[i, j, k];
                                        w[i, j, k] = w[i, j, k - 1];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


        /*****************************************************************************
         * Perform a single time step
         *****************************************************************************/
        public void time_step()
        {
            diffuse(u_old, ref u);
            diffuse(v_old, ref v);
            diffuse(w_old, ref w);

            project();

            Array.Copy(u, 0, u_old, 0, u.Length);
            Array.Copy(v, 0, v_old, 0, v.Length);
            Array.Copy(w, 0, w_old, 0, w.Length);

            advect(ref u, u_old, u_old, v_old, w_old, 2);
            advect(ref v, v_old, u_old, v_old, w_old, 3);
            advect(ref w, w_old, u_old, v_old, w_old, 4);

            project();

            Array.Copy(u, 0, u_old, 0, u.Length);
            Array.Copy(v, 0, v_old, 0, v.Length);
            Array.Copy(w, 0, w_old, 0, w.Length);
        }

        /*****************************************************************************
         * Solves the banded system given by the finite difference method applied
         * to the Poisson or diffusion equation using the iterative Gauss-Seidel method.
         * @inputs
         * double a - coefficient along diagonal entry
         * double[] c[6] - coefficient of all other nonzero entries in order from left 
         * to right (-k, -j, -i, +i, +j, +k)
         * double[, ,] b - right hand side
         * double[, ,] x0 - initial guess
         * out double[, ,] x1 - solution
         * 
         * TO DO: add Jacobi solver, which can be run on multiple cores
         ****************************************************************************/
        void gs_solve(double a, double[] c, double[, ,] b, double[, ,] x0, ref double[, ,] x1)
        {
            int Sx = x0.GetLength(0);
            int Sy = x0.GetLength(1);
            int Sz = x0.GetLength(2);

            //x1 = new double[Sx, Sy, Sz];

            int iter = 0;
            double res = 2 * solver_prams.tol;

            apply_boundary_conditions();

            while (iter < solver_prams.min_iter || (iter < solver_prams.max_iter && res > solver_prams.tol))
            {
                for (int k = 1; k < Sz - 1; k++)
                {
                    for (int j = 1; j < Sy - 1; j++)
                    {
                        for (int i = 1; i < Sx - 1; i++)
                        {
                            if (omega.obstacle_cells[i, j, k] == 0) //if cell not an obstacle
                            {
                                x1[i, j, k] = (b[i, j, k] - (c[0] * x0[i, j, k - 1] + c[1] * x0[i, j - 1, k] + c[2] * x0[i - 1, j, k] +
                                            c[3] * x0[i + 1, j, k] + c[4] * x0[i, j + 1, k] + c[5] * x0[i, j, k + 1])) / a;
                            }
                        }
                    }
                }

                apply_boundary_conditions();

                res = Utilities.compute_L2_difference(x0, x1);
                iter++;

                Array.Copy(x1, 0, x0, 0, x1.Length);
            }

            if (solver_prams.verbose)
            {
                Console.WriteLine("Gauss-Seidel solver completed with residual of {0} in {1} iterations", res, iter);
            }
        }
    }
}
