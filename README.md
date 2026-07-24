# SatEvol
## Satellite Array Evolution and Close Approach Tracking

SatEvol is a program that propagates a satellite population in low-Earth orbit (LEO) while tracking conjunctions. Input is provided as a list of Three-Line Element sets (TLEs) from which orbital elements are extracted and used to construct a Keplerian orbit for each satellite, with perturbations from Earth's oblateness taken into account (the J2 zonal coefficient). At each time-step, any objects closer than some threshold distance are logged and saved to a data file (format described below).

Plotting scripts are included to visualize the results, including conjunction tracks, relative distance and velocity histograms, time between conjunctions as a function of altitude, and a count of conjunctions over time. The simulation differentiates between satellites, debris, and derelict rocket bodies.

This code was produced primarily to verify close-encounter rates with the more heuristic [CRASH Clock](https://outerspaceinstitute.ca/crashclock/), a metric intended to provide the average time to first collision in the absence of any stationkeeping or collision avoidance maneuvers.

## Usage

There are a handful of parameters/flags at the beginning of `satEvol.py` to be configured:

`TLES_FILENAME [str]`:    Path to file containing list of TLEs to be read in, **or** or a checkpoint file from a previous simulation run (described below)

`OUT_SUFFIX [str]`:       String appended to all outfiles for a given run

`JD [float]`:             Julian date to propagate all satellites to before beginning the main integration. If any TLEs are provided with a date that differs from `JD` by more than `JDOFFSET [float]`, they will be omitted from the calculation.

`DT [float]`:             Length of time-step [s] (Satellites in LEO travel around 7-10 km/s, so a small timestep is important for cataloguing very close encounters)

`TTIME [float]`:          Total time of sim run [s]

`DCLOSE_METRE [float]`:   Minimum distance to catalogue a close encounter [m]

`READ_FROM_CHECKPOINT [bool]`:  Whether or not to read from a checkpoint file (True) or a list of TLEs (False)

`LINK [bool]`:  If true and reading from a list of TLEs, the code will attempt to find objects initialized within a small distance of each other (~100 m) and consolidate them into a single object. This is intended to handle modular spacecraft with multiple components that have different TLEs (e.g. the ISS), and uninitialized TLEs (which usually initialize at the origin).

`RANDOM_ORBITS [bool]`: If true, the RAAN, mean anomaly, and argument of periapsis (at epoch) for each object is randomized according to a uniform distribution between 0 and 2π before propagating. This is intended to simulate randomization of the orbital environment on long time-scales.

---

Once configured, the code can be run with

> `python3 satEvol.py`

While running the program displays the current simulation time and the number of satellites that currently fall within `DCLOSE_METRE` distance of eachother. Checkpoint files are periodically saved so that the simulation can be stopped and restarted. To do so, set `READ_FROM_CHEACKPOINT = True` and change the input filename to the most recent checkpoint file.

## Output

By default the outfile name is `track_out_` + OUT_SUFFIX + `.dat`, written to the `out` directory. The columns are:

`Time [s], Conjunction distance [m], Relative velocity [m/s], Altitude [m], Index of 1st satellite, Index of 2nd satellite, Name of 1st satellite, Name of 2nd satellite`

The helper script `plot_tracks.py` can be used to visualize the results. Simulation parameters need to be copied over to the relevant global variables before running. The script produces a plot of:

- Fitted conjunction distances over time ("tracks") for all orbiting objects
- 2D histogram of the minimum close approach distances and relative velocities for all objects
- 1D histogram of the relative velocities between all conjuncting objects
- 2D histogram of the time between conjunctions and the altitudes of all conjuncting objects
- Time series of conjunction frequency between the different object populations
- Frequency of close encounters as a function of altitude

## Dependencies

- NumPy v2.26+
- Pandas v2.3.1+
- SciPy v1.15.3+
- Numba v0.61.2 (Numba is still in active development, higher versions may or may not introduce issues)
- Python implementation of [SGP4 v2.24+](https://pypi.org/project/sgp4/)
- Matplotlib (for visualizing output)

`KeplerTools.py` is a helper module containing functions for common astrodynamical calculations. This is where the bulk of the numerical orbital propagation takes place.

`satArr.py` contains the class holding the collection of intialized satellites and their orbital elements. The class contains methods for propagating the orbits and figuring out which objects are in conjunction.

## Contact and Citation

For questions about the code, feel free to reach out to `skyeh@phas.ubc.ca`. 
If you use any of this code in your own work, or extend it, please cite [this paper](https://www.sciencedirect.com/science/article/pii/S0094576526004091).


