# Numerical Foundations

This chapter introduces three core classes of numerical methods that
appear throughout chemistry and many other quantitative applications:
root finding, optimisation, and the numerical solution of ordinary
differential equations. These techniques allow us to solve implicit
equations, determine unknown model parameters, locate minima on
potential-energy surfaces, and describe the time evolution of chemical
systems. 
Along the way, the chapter also introduces
central programming patterns for numerical work, including reusable
functions, iterative algorithms, and standard tools from `scipy`.

| Section | Covered Examples | New Concepts and Tools |
| --- | --- | --- |
| Root Finding | molar volumes from the Van der Waals equation | `while`-loops, default arguments, iterative methods, convergence criteria, finite-difference derivatives |
| Computational Optimisation | yield optimisation for a synthetic reaction | `for`-loops, vector-valued parameters, objective functions, termination vs. convergence, `scipy.optimize.minimize` |
| ODE Solvers | concentration profile modelling for a chemical reaction | 2D arrays, initial value problems, time-series plotting, `scipy.integrate.solve_ivp` |