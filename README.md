# `mdx-lib`: crystallographic computing in MATLAB.

**mdx-lib** is a general-purpose crystallographic computing library designed for processing and simulating _macromolecular diffuse X-ray scattering_ (MDX).

To see `mdx-lib` in action, check out [mdx-examples](https://github.com/ando-lab/mdx-examples).

## Manuscripts

[Meisburger, Case, & Ando, 2020]: https://doi.org/10.1038/s41467-020-14933-6
Meisburger, SP, Case, DA & Ando, N. (2020). Diffuse X-ray scattering from correlated motions in a protein crystal. _Nat Commun_ **11**, 1271. <https://doi.org/10.1038/s41467-020-14933-6>

[Meisburger, Case, & Ando, 2022]: https://doi.org/10.1101/2022.08.22.504832
Meisburger, SP, Case, DA & Ando, N. (2022). Robust total X-ray scattering workflow to study correlated motion of proteins in crystals. _bioRxiv_ 2022.08.22.504832; doi: <https://doi.org/10.1101/2022.08.22.504832>


## Version History

### Version 1.2

Added scripts to support GOODVIBES and DISCOBALL analysis described in [Meisburger, Case, & Ando, 2022].

- New tools to represent atomic scattering factors and ADPs `latt.Blob`
- New elastic network modelling tools `nm.ElasticNetwork_v2`, `nm.Cell_v2`, and `nm.LatticeDynamics`. These will replace old versions in future release.
- Demos moved to separate repository [mdx-examples](https://github.com/ando-lab/mdx-examples).
- High-level scripts to process maps `proc.script.MapTools`, `proc.script.GridDesigner`, `proc.script.MapInterp`
- Scripts for GOODVIBES analysis: `proc.script.ElasticNetworkTools`, `proc.script.LatticeDynamicsTools`
- Scripts for DISCOBALL analysis `proc.script.DeltaPDFTools`
- Scripts for wrangling crystallographic data files: `proc.script.ImportMTZ`, `proc.script.ImportPDB`

### Version 1.1

Released at the same time as [Meisburger, Case, & Ando, 2020] in a separate branch (natcomm). Includes scripts to reproduce the data processing and model fitting of lysozyme (P1 space group). There were also extra atomistic modeling and normal mode analysis libraries not included in the master branch.

### Version 1.0

Released at the same time as [Meisburger, Case, & Ando, 2020]. This more stable version includes all of the libraries for data processing, but not for atomistic modeling or normal mode analysis.

####

## Repository Contents

- `+geom` - diffraction geometry and scattering corrections
- `+grid` - fractional Miller index grids for voxel integration
- `+io` - read and write files
- `+latt` - direct and reciprocal lattices, Fourier transforms, electron density of atoms, atomic displacement parameters
- `+model` - atomistic modeling, material properties, coherent and incoherent scattering factors
- `+nm` - crystalline elastic network modeling and normal mode calculations
- `+proc` - algorithms, high-level scripts, and batch processing routines
- `+symm` - space groups, symmetry operators, and affine transformations
- `+util` - functions and class definitions used by other libraries

## Requirements

### Hardware

The data processing libraries are designed to run on a desktop computer with a large amount of RAM and multiple CPU cores for parallel processing. RAM is needed mostly for processing large datasets or if parallel Batch processing is used. It is recommended to run small-scale tests before scaling up.

The development system was a Mac with 64 GB RAM and 10 cores @ 3 GHz/core.

### Software

`mdx-lib` runs within MATLAB. For full functionality, the following toolboxes should be installed:

 - Parallel Processing Toolbox (The Mathworks)
 - Statistics and Machine Learning Toolbox (The Mathworks)

The software has been tested on MATLAB version R2021a running on macOS Version 10.14.6 (Mojave). Unix-like paths are assumed throughout (Windows is not tested or supported currently).

## Installation Instructions

Download a release or clone this repository.

Add mdx-lib directory to the MATLAB path. This can be done using the path editor, or for each session by typing (command prompt):

```matlab
>>> addpath('/<path>/<to>/mdx-lib')
```

### Installation issues

#### cbf decompression fails

CBF decompression is accelerated using MATLAB coder to translate the m-file [cbf_decompress.m](+io/@DectrisCBF/private/cbf_decompress.m) into c and compile it into a mex file:
```
+io/@DectrisCBF/private/cbf_decompress_mex.mexmaci64
```
Since the mex file was compiled on an intel-based Mac, it will only work on similar systems. If using a different OS or architecture, cbf_decompress will need to be re-compiled. The project file is included here:
```
+io/@DectrisCBF/private/cbf_decompress.prj
```

Users have also reported the following error:

> cbf_decompress_mex.mexmaci64 cannot be opened because the developer cannot be verified

Running this in the terminal should fix the issue (substitute with path to repository on your machine):

```bash
sudo xattr -r -d com.apple.quarantine /<path>/<to>/mdx-lib

sudo find /<path>/<to>/mdx-lib -name \*.mexmaci64 -exec spctl --add {} \;
```

## mdx2 ...

We're working on a complete re-implementation in python! Check it out: [mdx2](https://github.com/ando-lab/mdx2)
