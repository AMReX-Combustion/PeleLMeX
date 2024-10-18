## EB\_ODEQty
A 2D inflow/outflow setup with an optional EB cylinder in the middle of the flow. Demonstrates how to use ProblemSpecificFunctions and the ODE quantities.  The ODE quantities experience simple exponential decay that gets stiffer for each quantity.  Specifically, $\frac{\partial B_k}{\partial t} = \gamma \cdot 10^{k+1} B_k$, for $k = 0, 1,\dots,$ NUM_ODE and $\gamma < 0$.
