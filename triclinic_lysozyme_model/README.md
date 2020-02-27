# Triclinic Lysozyme Model

Create vibrational models of triclinic lysozyme and fit them to measured diffuse scattering and atomic displacement parameters (ADPs).

## Input data

The scripts require `mdx-lib` and helper functions in `model/code` and `fit/code` as well as the following input data files:

* Refined atomic coordinates and ADPs (`.pdb`)
* Measured structure factors (`.mtz`)
* Diffuse scattering map (`.h5`)

For convenience, the atomic coordinates ([PDB ID 6o2h](https://www.rcsb.org/structure/6O2H)) are automatically downloaded to `model/6o2h.pdb` when running `get_atomic_model` (see below). Measured structure factors are included in `model/6o2h_aimless_truncate.mtz` (diffraction images were integrated with [xds](http://xds.mpimf-heidelberg.mpg.de/) and scaled/merged with aimless and truncate in [ccp4](http://www.ccp4.ac.uk/)). The diffuse scattering map can be obtained from the raw diffraction images (see [`triclinic_lysozyme_maps`](../triclinic_lysozyme_maps)) or downloaded from the CXIDB ([ID 128](https://www.cxidb.org/id-128.html)). In the latter case, `setup_environment.m` must be edited with the path to the h5 file.

## Instructions

First, inspect `setup_environment.m` to verify the paths to input data and code are correct for your system. Note that the contents of this folder can be moved to a different location before running the scripts, as long as the paths in `setup_environment.m` are modified accordingly.

Next, the scripts in the main directory should be inspected and run in the following order:

1. `setup_environment`
1. `get_atomic_model`
1. `elastic_network_model`
1. `lattice_dynamics_model`
1. `reference_halos`
1. `reference_calc`
1. `fit_lattice_dynamics_model_to_reference`
1. `fit_elastic_network_model_to_adps`

Each script saves its results to a `.mat` of the same name in the directories `model` and `fit`.
