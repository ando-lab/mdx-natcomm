# Triclinic Lysozyme Maps

Create maps of diffuse scattering from triclinic lysozyme on an absolute intensity scale

## Input data

The scripts require `mdx-lib` and helper functions in `calc/code` and `export/code` as well as the following input data:

* Diffraction images (`.cbf`)
* Diffraction geometry from [xds](http://xds.mpimf-heidelberg.mpg.de/) (`INTEGRATE.LP` and `XDS.INP`)
* Measured Bragg peak intensities (`.mtz`)

The diffraction images can be downloaded from the SBGrid Data Bank ([ID 747](https://data.sbgrid.org/dataset/747/)) to any convenient location. To obtain the diffraction geometry and Bragg peak intensities, these images were integrated with [xds](http://xds.mpimf-heidelberg.mpg.de/) and scaled/merged with aimless in [ccp4](http://www.ccp4.ac.uk/). It is not necessary to re-run XDS, as the required `.INP`, `.LP` and `.mtz` files are included in within the `proc` directory.

## Instructions

First, inspect `setup_environment.m` to verify the paths to input data and code are correct for your system. Note that the contents of this folder can be moved to a different location before running the scripts, as long as the paths in `setup_environment.m` are modified accordingly.

Next, the scripts in the main directory should be inspected and run in the following order:

1. `setup_environment`
1. `process`
1. `calc_half_integer_map`
1. `calc_intensity_statistics`
1. `calc_variational_intensity`
1. `calc_interpolated_7x7x7`
1. `export_maps`

Data reduction occurs in the `process` script, which produces output `.mat` and `.log` files in the `proc` directory. It is designed to be run step-wise in cell mode in the MATLAB editor. The main result is the fine diffuse map, which is stored in `proc/mdx/mergeFine.mat`.

The `calc_*` scripts further process this map, each producing an output `.mat` file of the same name in the `calc` directory. Finally, `export_maps` writes the various maps to a single HDF5 file `export/triclinic_lysozyme_maps.h5`.  This output file can be compared with the author's version deposited in the CXIDB ([ID 128](https://www.cxidb.org/id-128.html)).
