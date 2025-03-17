
# SatEvol
## Satellite Cloud Evolution and Close Approach Tracking

This code reads satellites from a list of TLEs and propagates them in orbit around the Earth using the SGP4 model. At each time-step all satellites are checked for close approaches with other satellites and logged.

## To-do List:

In approximately decreasing priority,

1. ~~Implement distance more efficient nearest-neighbour search - SciPy KDTree should be O(n log n) instead of O(n^2).~~
2. Vectorize KeplerTools
> Partially done
3. Implement numba (JIT compilation) and (probably) parallelize
> Also partially done, only in KT right now, will need to reorganize program structure to parallelize better
4. Extend SatEvol to examine phase-space mixing.
4. Re-evaluate performance and decide if it needs to be translated
5. Rewrite all I/O with list comprehensions (or array comprehensions if it helps significantly)
6. Remove plotting and create a helper program to take care of it afterwards for deployment

## Plots and Quantities Needed:

1. Phase-space mixing histogram for starlink constellation. Sampled ~hourly. 2D histogram?
2. Time scale for close encounters as a function of altitude.
3. (Maybe) Object threshold for avoidance maneuver cascade

## Changelog:

- 11/03/25: Implemented KDTree distance search
- 12/03/25: Optimized a chunk of KeplerTools, implemented numba for jit compilation.
- 17/03/25: Minor file I/O clean-up, including outfiles for 4-orbit test run.