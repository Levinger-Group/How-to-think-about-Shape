# Coordinate-Pair Eccentricity (CPE)

Written by CG, 2021

Code needed to compute the CPE. Calculating the CPE itself is trivial but computing the semi-axes is slightly more computationally expensive. This code assumes you are working with a molecular dynamics trajectory. Options for creating semi-axes are to read in an xvg file from a GROMACS gyrate analysis, compute the semi-axis of a single frame of the trajectory, compute the semi-axis of all frames in the trajectory in serial, and compute the semi-axis of all frames in the trajectory in parallel. 

Note that multiprocessing (via Dask) is only available for mdanalysis >v2.0
