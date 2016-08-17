Title: FFD library readme file
Author: Lukas Bystricky
Date: August 17th, 2016

# FAST FLUID DYNAMICS LIBRARY

This is a C# implementation of the fast fluid dynamics algorithm (FFD) described by [Jos Stam](http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/ns.pdf). This code has been written as part of the Google Summer of Code with Christoph Waibel of [EMPA](https://www.empa.ch/web/empa/) as a mentor. The goal of this project was to use FFD to quickly analyze the effects of urban windflow on natural ventilation of buildings. 

## Algorithm outline

FFD solves the incompressible Navier-Stokes equations using a fully implicit projection method. The incompressible Navier-Stokes equations are given by:

![Navier-Stokes](/img/navier-stokes.png)

There are 3 main steps involved in FFD. Starting at $\mathbf{u}^n approx \mathbf{u}(t_0)$ we solve for $\mathbf{u}^{n+1} \approx \mathbf{u}(t_0 + \delta)$ as follows:

1. Resolve the diffusion term using implicit Euler:

$$\frac{\mathbf{u}^* - \mathbf{u}^n}{\delta} = \nu\Delta\mathbf{u}^*$$