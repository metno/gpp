Gridded post-processor
======================

.. image:: https://travis-ci.org/metno/gridpp.svg?branch=master
    :target: https://travis-ci.org/metno/gridpp

.. image:: https://coveralls.io/repos/metno/gridpp/badge.svg?branch=master&service=github
    :target: https://coveralls.io/github/metno/gridpp?branch=master 

Gridpp is a command-line tool that post-processes weather forecasts in NetCDF format. The program
performs two types of post-processing: Downscaling and calibration. Gridpp downscales forecast from
a coars grid to a finer grid using a variety of interpolation methods. Gridpp then calibrates the
forecasts by applying corrections to each gridpoint. Gridpp is modular, so any combination of
downscaling and calibration can be selected.

For information on how to use the software, check out the wiki page:
https://github.com/metno/gridpp/wiki


Features
--------

* Downscaling from one grid to another using methods such as nearest neighbour, bilinear
  interpolation, correction based on elevation gradients, and more.
* Calibration such as bias-correction, neighbourhood smoothing, probabilistic calibration, and more.
* Support for Netcdf files with great flexibility in variable names, number of dimensions, and names
  of dimensions.

Installing on Ubuntu
---------------------

**Install dependencies**

.. code-block:: bash

  sudo apt-get update
  sudo apt-get install libboost-all-dev
  sudo apt-get install libgsl0-dev libblas-dev
  sudo apt-get install netcdf-bin libnetcdf-dev
  sudo apt-get install cmake
  sudo apt-get install libarmadillo6 libarmadillo-dev

Optional: If you are developing the code, install the google test library:

.. code-block:: bash

  sudo apt-get install libgtest-dev
  sudo apt-get install lcov
  sudo sudo gem install coveralls-lcov

**Installing from source**

1. Download the source code from the latest release: https://github.com/metno/gridpp/releases. Unzip
   the file and navigate into the extracted folder.

2. Create a build directory

.. code-block:: bash

   mkdir build

3. Run cmake to set up installation

.. code-block:: bash

   cd build
   cmake ..

This will set up gridpp to be installed in /usr/local/bin on Ubuntu. To specify a custom installation path, use:

.. code-block:: bash

   cmake .. -DCMAKE_INSTALL_PREFIX=<custom path>

If you are developing code, run the following instead:

.. code-block:: bash

   cmake .. -DCMAKE_BUILD_TYPE=DEBUG -DENABLE_TESTS=ON -DGTEST_DIR=/usr/src/gtest

where ``/usr/src/gtest`` is the location of the google test library code.

4. Run cmake to install

.. code-block:: bash

   cmake --build .
   cmake --build . --target install

**Installing from debian packages**

Follow instructions here: https://launchpad.net/~metno/+archive/ubuntu/gridpp


Copyright and license
---------------------
Copyright (C) 2014-2018 MET Norway. Gridded post-processor is licensed under `GPL
version 2 <https://github.com/metno/gridpp/blob/master/LICENSE>`_ or (at
your option) any later version.

Contact
-------
| E-mail: `thomasn@met.no <mailto:thomasn@met.no>`_
