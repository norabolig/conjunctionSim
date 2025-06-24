# SatEvol
## Satellite Array Evolution and Close Approach Tracking

SatEvol is a program that propagates a satellite population in low-Earth orbit (LEO) while tracking conjunctions. Input is provided as a list of Two/Three-Line Element sets (TLEs) from which orbital elements are extracted and used to construct a Keplerian orbit for each satellite, with perturbations from apsidal and nodal precession taken into account. At each time-step any conjuncting objects closer than a threshold distance are logged and saved to a data file (format described below).

Plotting scripts are included to visualize the results, including conjunction tracks, relative distance and velocity histograms, time between conjunctions as a function of altitude, and a count of conjunctions over time. The simulation differentiates between satellites, debris, and derelict rocket bodies.

## Usage

There are a handful of parameters at the beginning of `satEvol.py` to be configured:


`TLES_FILENAME [str]`:    Path to file containing list of TLEs to be read in

`PHASE_FOUT [str]`:       Name for outfile containing initial phase-space configuration (i.e. positions and velocities for every satellite at t=0)

`OUT_SUFFIX [str]`:       String appended to all outfiles for a given run

`DT [float]`:             Length of time-step [s]

`TTIME [float]`:          Total time of sim run [s]

`PHS_INT_S [float]`:      How long to wait between phase-space snapshots if `EXAMINE_PHS=True`

`DCLOSE_METRE [float]`:   Minimum distance to catalogue a close encounter [m]

Once configured, the code can be run with

> `python3 satEvol.py`

While running the program displays the current time-step and the number of satellites that currently fall within `DCLOSE_METRE` distance of eachother. 
By default the outfile name is `track_out_` + OUT_SUFFIX + `.dat`, written to the `out` directory. The columns are:

`[Time], [Conjunction distance], [Relative velocity], [Altitude], [Index of 1st satellite], [Index of 2nd satellite], [Name of 1st satellite], [Name of 2nd satellite]`

## Planned Features

- System to handle collisions and inject debris into the simulation.
- Wider range of perturbative effects (for long time-scale runs)
- Move number crunching to a compiled language or the GPU.

## To-do List:

In approximately decreasing priority,

- ~~Do run, count all instances of close approaches under 1km, 2km, 4km, 8km, etc... (up to 15km)~~
- ~~Implement mode to calculate how long until first <200m close approach~~

- Make plots of orbital elements for debris / full-catalogue / starlinks
- ~~Implement apsidal precession~~
- Implement check-pointing system

1. ~~Implement distance more efficient nearest-neighbour search - SciPy KDTree should be O(n log n) instead of O(n^2).~~
2. ~~Vectorize KeplerTools~~
3. ~~Implement numba (JIT compilation)~~ and (probably) parallelize
4. Re-evaluate performance and decide if it needs to be translated
5. Remove plotting and create a helper program to take care of it afterwards for deployment

## Changelog:

- 11/03/25: Implemented KDTree distance search.
- 12/03/25: Optimized a chunk of KeplerTools, implemented numba for jit compilation.
- 17/03/25: Minor file I/O clean-up, including outfiles for 4-orbit test run.
- 23/03/25: Implemented phase-space mixing. Fixed bug in outfiles.
- 26/03/25: Performed 48 hour test run, dt = 0.05s w/ 15km threshold.
- 02/04/25: Cleaned up some memory management.
- 18/06/25: Implemented apsidal precession
