
# SatEvol
## Satellite Cloud Evolution and Close Approach Tracking

This code reads satellites from a list of TLEs and propagates them in orbit around the Earth using the SGP4 model. At each time-step all satellites are checked for close approaches with other satellites and logged.

## To-do List:

In approximately decreasing priority,

1. Implement distance more efficient nearest-neighbour search.
2. Implement numba (JIT compilation) and parallelize
3. Rewrite all I/O with list comprehensions (or array comprehensions if it helps significantly)
4. Re-evaluate performance and decide if it needs to be translated

The code is almost entirely in "raw" Python right now using lists. SciPy and numpy should contain a significant portion of these operations
as precompiled nicely vectorized C/Fortran code; should save the headache of doing a full translation. The following might be useful: