# SatEvol
## Satellite Cloud Evolution and Close Approach Tracking

This code reads satellites from a list of TLEs and propagates them in orbit around the Earth using the SGP4 model. At each time-step all satellites are checked for close approaches with other satellites and logged.

## How to run:

There are a handful of parameters at the beginning of `satEvol.py` to be configured:


`TLES_FILENAME [str]`:    Path to file containing list of TLEs to be read in

`PHASE_FOUT [str]`:       Name for outfile containing initial phase-space configuration (i.e. positions and velocities for every satellite at t=0)

`OUT_SUFFIX [str]`:       String appended to all outfiles for a given run

`EXAMINE_PHS [bool]`:     Whether or not examine phase-space mixing with snapshots

`PLOT [bool]`:            Wether or not to create plots directly with satEvol.  

`DT [float]`:             Length of time-step [s]

`TTIME [float]`:          Total time of sim run [s]

`PHS_INT_S [float]`:      How long to wait between phase-space snapshots if `EXAMINE_PHS=True`

`DCLOSE_METRE [float]`:   Minimum distance to catalogue a close encounter [m]

Once configured, the code can be run with

> `python3 satEvol.py`

While running the program displays the current time-step and the number of satellites that currently fall within `DCLOSE_METRE` distance of eachother. 
By default the outfile name is `track_out_` + OUT_SUFFIX + `.dat`, written to the `out` directory. The columns are:

`[Time], [Conjunction distance], [Relative velocity], [Altitude], [Index of 1st satellite], [Index of 2nd satellite], [Name of 1st satellite], [Name of 2nd satellite]`

## To-do List:

In approximately decreasing priority,

- Do run, count all instances of close approaches under 1km, 2km, 4km, 8km, etc... (up to 15km)
- Implement mode to calculate how long until first <200m close approach

1. ~~Implement distance more efficient nearest-neighbour search - SciPy KDTree should be O(n log n) instead of O(n^2).~~
2. Vectorize KeplerTools
> Partially done
3. Implement numba (JIT compilation) and (probably) parallelize
> Also partially done, only in KT right now, will need to reorganize program structure to parallelize better
4. Extend SatEvol to examine phase-space mixing.
    - Write nearest-neighbour quantities to outfile
    - Constrain distances under X km to examine tail (handled by plotting helper?)
    - Add counter for time
    - Make histograms 2d, add lines on 1d to indicate median quantities

4. Re-evaluate performance and decide if it needs to be translated
5. Remove plotting and create a helper program to take care of it afterwards for deployment

## Plots and Quantities Needed:

1. Phase-space mixing histogram for starlink constellation. Sampled ~hourly. 2D histogram?
2. Time scale for close encounters as a function of altitude.
3. (Maybe) Object threshold for avoidance maneuver cascade

## Changelog:

- 11/03/25: Implemented KDTree distance search.
- 12/03/25: Optimized a chunk of KeplerTools, implemented numba for jit compilation.
- 17/03/25: Minor file I/O clean-up, including outfiles for 4-orbit test run.
- 23/03/25: Implemented phase-space mixing. Fixed bug in outfiles.
- 26/03/25: Including outfiles for 48 hour test run, dt = 0.05s w/ 15km threshold.