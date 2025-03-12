
# SatEvol
## Satellite Cloud Evolution and Close Approach Tracking

This code reads satellites from a list of TLEs and propagates them in orbit around the Earth using the SGP4 model. At each time-step all satellites are checked for close approaches with other satellites and logged.

## To-do List:

In approximately decreasing priority,

1. ~~Implement distance more efficient nearest-neighbour search - SciPy KDTree should be O(n log n) instead of O(n^2).~~
2. Vectorize KeplerTools
3. Implement numba (JIT compilation) and (probably) parallelize
4. Re-evaluate performance and decide if it needs to be translated
5. Rewrite all I/O with list comprehensions (or array comprehensions if it helps significantly)
6. Remove plotting and create a helper program to take care of it afterwards for deployment

## Other Notes:

- Nothing rn :3