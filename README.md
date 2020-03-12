# gr_calc_matlab
A simple symbolic general relativity calculator for Matlab.

Instantiate a class with "gr_calc(metric,coordinates)" and use the various functions such as riemann, christoffel, ricci, einstein.

See the live script for an example.

## Simulink blocks

A simulink solver for geodesic equations. The metrics can be generated with matlabFunction.

The solver is for the schwarchild metric.  Other metrics can be solved by 

1.) generating the inverse metric and the inverse derivative metric

2.) adding appropriate selectors and/or signal seperators into and attaching the right signals to the metric

3.) add any necessary constant blocks


Also, be careful with initial conditions.