# Gridded post-processor

[!["Latest release"](https://img.shields.io/github/v/release/metno/gridpp.svg)](https://github.com/metno/gridpp/releases)
[![C/C++ CI](https://github.com/metno/gridpp/workflows/C/C++%20CI/badge.svg)](https://github.com/metno/gridpp/actions)

Gridpp a is post-processing tool for gridded weather forecasts. It consists of a **library** of commonly-used methods and a **command-line tool** that applies these methods to forecast fields in NetCDF files.

Gridpp is written in C++ but offers python bindings to the functions in the library. The tool is used at MET Norway to produce operational weather forecasts for Yr (https://www.yr.no).

Gridpp is currently under active development and the current version is a prototype for testing. Feedback
is welcome, either by using the issue tracker in Github, or by contacting Thomas Nipen (thomasn@met.no).

## Documentation
For information on how to use gridpp, check out the wiki at https://github.com/metno/gridpp/wiki.

## Features
- Methods for **downscaling** a forecast from a coarse grid to a fine grid
- Methods for **calibrating** a downscaled grid, such as quantile mapping
- Computationally efficient **neighbourhood** methods to compute neighbourhood min, mean, max, and any quantile.
- Data assimilation using **optimal interpolation** (OI) to merge observations and gridded forecasts (deterministic or ensemble)
- Efficient data structures for nearest location lookup in a vector or grid of locations
- Command-line client with support for Netcdf files with flexibility in how variables and dimensions are configured

## Example

The following computes a moving neighbourhood mean with a half-width of 7 gridpoints (i.e. 15x15 neighbourhood)

```python
import gridpp
import numpy as np

field = np.random.rand(300, 200)
halfwidth = 7
gridpp.neighbourhood(field, halfwidth, gridpp.Mean)
```

![Example](extras/image.gif)

## Required dependencies
- [Boost](https://www.boost.org/) >= 1.59
- [Armadillo](http://arma.sourceforge.net/) >= 6.6
- [GNU Scientific Library](https://www.gnu.org/software/gsl/)
- [Netcdf](https://www.unidata.ucar.edu/software/netcdf/)

On Ubuntu Bionic, these can be installed like this:
```bash
sudo apt-get update
sudo apt-get install libboost-all-dev
sudo apt-get install libgsl0-dev libblas-dev
sudo apt-get install netcdf-bin libnetcdf-dev
sudo apt-get install libarmadillo6 libarmadillo-dev
```

Note that Ubuntu Xenial only has Armadillo 6.5 in its apt repository. In that case you need to install  [Armadillo 6.6](http://arma.sourceforge.net/) or later manually.

## Installing the python bindings from pip

The easiest is to install the latest release of the package using pip. Provided you have installed the dependencies listed above, you can install the most recent release of the python package as follows:
```bash
pip3 install gridpp --user
```

To check that the installation worked, run the following in python3:
```python
import gridpp
print(gridpp.version())
```

## Full gridpp installation from source

1. Either download the source code from the [latest release](https://github.com/metno/gridpp/releases), unzip
   the file and navigate into the extracted folder; or clone the repo from github.

2. Install extra requirements

These are only required when installing from source
```
sudo apt install swig cmake
```

3. Set up cmake installation

```bash
mkdir build
cd build
cmake ..
```

4. Install the C++ library

```bash
sudo make install
```
This will install the library in `/usr/local/lib/libgridpp.so` and the gridpp command-line client in
`/usr/local/bin/gridpp`. To specify a custom installation path, use the following in step 3:

```bash
cmake .. -DCMAKE_INSTALL_PREFIX=<custom path>
```

5. Install the python bindings

```bash
make install-python-user
```

This installs the python bindings in
`~/local/lib/python3.6/site-packages/gridpp.py`. To install the python bindings system-wide, use `sudo make install-python` instead.

6. Install the R bindings

```bash
make build-r
```

Currently, the R package is not installed centrally, but instead is placed in `extras/SWIG/R/gridpp.R` in the build directory.

## Copyright and license
Copyright © 2014-2020 Norwegian Meteorological Institute. Gridpp is licensed under the GNU LEsser General
Public License (LGPL). See LICENSE file.

## Contact
E-mail: Thomas Nipen (thomasn@met.no)
