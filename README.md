# FAST FLUID DYNAMICS LIBRARY

This is a C# implementation of the fast fluid dynamics algorithm (FFD) described by [Jos Stam](http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/ns.pdf). This code has been written as part of the Google Summer of Code with Christoph Waibel of [EMPA](https://www.empa.ch/web/empa/) as a mentor. The goal of this project was to use FFD to quickly analyze in 3D the effects of urban windflow on natural ventilation of buildings. 


## How to Use

In the Test folder there are 3 drivers that demonstrate how to set up various simulations:

* `CavityDriver.cs` simulates a lid driven cavity, a well known CFD benchmark problem
* `ConvergenceDriver.cs` runs a convergence study using the exact 3D unsteady solutions given by [Ethier and Steinman](http://onlinelibrary.wiley.com/doi/10.1002/fld.1650190502/abstract)
* `BuildingSimulation.cs` runs a simulation of windflow past 3 buildings

![Lid driven cavity](img/cavity_re100_mag.png)![Lid driven cavity](img/cavity_re100_u.png)

Cross section of velocity magnitude and x component of velocity for lid driven cavity at RE=100. 

![Streamlines past buildings](img/flow_multiple_buildings.png)
Streamlines of wind flow past 3 buildings.

## Notes on Implementation

In Stam's original implementation (and in numerous others available online) cell centred finite difference is used to discretize the equations. In this implementation however we use staggered grid finite difference, which is the standard finite difference implementation in CFD. This is done to prevent spurious pressure oscillations near the boundary which can occur in cell centred finite diference for the Navier-Stokes equations. This does not change much in the algorithm or solvers, but makes enforcing the boundary conditions significantly more complicated. 

The only linear solver implemented is a simple Jacobi solver, as was the case in Stam's original implementation. This was done for ease of programming and the potential for easy parallelization, but should be replaced for any serious engineering studies as it converges far too slowly to be useful (see "future work"). 

## Code Description

### `FluidSolver.cs`

The FluidSolver class contains all the functions needed to advance a simulation 1 time step. To initialize it you must provide:

1. The initial velocity in each coordinate direction (`u0[, ,]`, `v0[, ,]` and `w0[, ,]`)
2. A domain (inherited from `Domain.cs`) that contains information about the mesh, obstacles and boundary conditions
3. A time step size
4. A viscosity
5. A `solver_struct` that contains:
    * solver tolerance
    * maximum and minimum number of iterations
    * backtrace order (1 or 2)
    * flag to indicate verbose console output

Since we are using a staggered grid for our finite difference scheme the sizes of the arrays can be confusing. Let us consider a domain which we wish to split into `Nx` cells in the x-direction, `Ny` cells in the y-direction and `Nz` cells in the z direction. Since we have a layer of ghost cells around our domain, in order to construct such a domain we would need a computational domain with `Nx+2`, `Ny+2` and `Nz+2` cells in the x, y and z directions. The centre of each of these cells will contain a pressure value.

The velocity values are defined on middle of the cell faces normal to their direction. This means that our velocity arrays will have the following sizes:

* `u`: `Nx+1`, `Ny+2`, `Nz+2`
* `v`: `Nx+2`, `Ny+1`, `Nz+2`
* `w`: `Nx+2`, `Ny+2`, `Nz+1`

### `Domain.cs`

The Domain class contains information about the computaional domain such as the domain size, the mesh size, obstacles, boundary conditions and exact solutions (if they exist). As mentioned above, in order to properly handle boundary conditions, there is a layer of ghost cells around the domain. These cells are marked as obstacle cells, even if the boundary is not a physical obstacle (i.e. inflow or outflow). The function `set_ghost_flags()` sets these cells and flags adjacent cells as the appropriate boundary cells when the Domain is constructed. 

The function `add_obstacle(xmin, xmax, ymin, ymax, zmin, zmax)` marks all cells inside the box [`xmin`, `xmax`] X [`ymin`, `ymax`] X [`zmin`, `zmax`] as obstacles and flags adjacent cells as the the appropriate boundary.

If an exact solution is known, the function `exact_solution(...)` in the base class can be overwritten to compute it. This can be used to perform convergence studies.

### `PostProcessor.cs`

The PostProcessor class contains routines to export the data (velocity, pressure, errors, geometry information) to [VTK](http://www.vtk.org/) files. This data can then be analyzed or plotted using such programs as [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/) or [ParaView] (http://www.paraview.org/). 

### `Utlilities.cs`

The Utilities class is a static class that contains general routines that may be used by multiple classes. These are:
* trilinear interpolation
* L2 difference between vectors
* L2 and L-infinity error calculation
* check if point inside obstacle/ghost cell

## Future Work

### Improved linear solvers

In many existing implementations of FFD, simple iterative solvers like Jacobi, Gauss-Seidel or Successive over-relaxation were used. This implementation uses a Jacobi solver because of the fact that it is easy to implement and easy to parallelize. Unfortunately this solver converges very slowly. This may not have been an issue for video game applications since accuracy isn't necessarily important, but for engineering applications on large domains this can be a big problem. A conjugate gradient or multigrid solver would perform much better. 

### Boussinesq approximation

In urban windflow, buoyancy effects are known to be very important. An extension of this model could implement the Boussinesq approximation to account for buoyancy. 

