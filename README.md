# `mdx-lib`

## Description

`mdx-lib` is a MATLAB library for processing and simulating macromolecular diffuse X-ray scattering (MDX). See our paper [1] for details.

## Repository Contents

- `+io` - read and write files (e.g. from XDS, aimless, Pilatus detectors)
- `+geom` - diffraction geometry and scattering corrections
- `+grid` - fractional Miller index grids for voxel integration
- `+proc` - batch processing scripts for primary data reduction
- `+latt` - direct and reciprocal lattices, Fourier transforms
- `+symm` - space group symmetries
- `+util` - functions and class definitions used by other libraries
- `+model` - atomic coordinates and scattering factors
- `+nm` - rigid-body crystalline elastic network for scattering simulation
- `demo` - a small simulated dataset and MATLAB scripts to process it. See Demonstration, below.
- `triclinic_lysozyme_maps` - scripts to process diffuse scattering from triclinic lysozyme data as described in Ref. [1]
- `triclinic_lysozyme_model` - scripts that fit elastic network and lattice dynamics models of triclinic lysozyme as described in Ref. [1]

## System Requirements

### Hardware Requirements

The processing libraries are designed to run on a desktop computer with a large amount of RAM and multiple CPU cores for parallel processing. The development system had the following specs:

- RAM: 64 GB
- CPU: 10 cores, 3 GHz/core

The demo uses simulated data with much small images, and will run quickly on a modest machine without parallel processing (see below).

### Software dependencies

`mdx-lib` runs within MATLAB. For full functionality, the following toolboxes should be installed:

 - Parallel Processing Toolbox (The Mathworks)
 - Statistics and Machine Learning Toolbox (The Mathworks)

The software has been tested on MATLAB version R2018a running on macOS Version 10.14.6 (Mojave).

## Installation Guide

Download the repository and add its directory to the MATLAB path.

## Demonstration

We have prepared a small demonstration in the `demo/` folder, which contains the following:

- `setup_environment.m` - function to set up MATLAB paths
- `process.m` - integrate, scale, and merge the simulated dataset
- `expected_results/` - published output of `process.m` in html format.
- `images/` - simulated diffraction images
- `code/` - helper functions

The demonstration can be run as follows:
1. Navigate to the `demo/` folder in MATLAB
2. Run `setup_environment` to set the paths
3. Run `process.m`. The script is designed to run in cell mode (MATLAB editor), but it can also be executed directly from the command line.

The demo creates a subdirectory `proc/`, where intermediate results and log files are stored, and produces several figures (see `expected_results/`).

The run time on a normal desktop computer is about 5 minutes.

For further information, see comments in `process.m`.

## References

[1]: Meisburger SP, Case DA, Ando N. (2020). Diffuse X-ray Scattering from Correlated Motions in a Protein Crystal. Nature Communications. http://dx.doi.org/10.1038/s41467-020-14933-6
