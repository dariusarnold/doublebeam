# Double gaussian beam

Find out fracture orientation and spacing from seismic data.

## Setup

Installation on Debian 10:

- Install g++, gcc, cmake and git
- Install libboost-program-options1.67-dev, libeigen3-dev, libfmt-dev, libfftw3-dev, libmsgsl-dev
- Clone the repository at https://github.com/dariusarnold/double-beam
- Change into double-beam folder
- Create a build directory and change to it
- Call cmake .. to create the Makefile and check if all required libraries are installed
- If you get the error 
  CMake Error at /usr/share/cmake-3.13/Modules/FindPackageHandleStandardArgs.cmake:137 (message):
  Could NOT find ZLIB (missing: ZLIB_LIBRARY ZLIB_INCLUDE_DIR)
  Call Stack (most recent call first):
  /usr/share/cmake-3.13/Modules/FindPackageHandleStandardArgs.cmake:378 (_FPHSA_FAILURE_MESSAGE)
  /usr/share/cmake-3.13/Modules/FindZLIB.cmake:114 (FIND_PACKAGE_HANDLE_STANDARD_ARGS)
  external/cnpy/CMakeLists.txt:12 (find_package)
  Install zlib1g-dev
- Call make doublebeam to create the doublebeam executable
- For info about how to run the doublebeam program, consult ./doublebeam -h

Required libraries and their Debian package names if applicable.

- fmt (libfmt-dev)
- Eigen3 (libeigen3-dev)
- FFTW3 (libfftw3-dev)
- Microsoft Guideline Support Library (libmsgsl-dev)
- Boost (libboost-dev)

OpenMP is optional. If it is found, parallel processing will be enabled. 

## Configuration

There are two configuration settings: ``config.cpp`` and the config file containing parameters for 
the algorithm. ``config.cpp`` provides customisation towards the general setup such as path names 
and typical filenames. Typically you would only have to change it once to adapt to your naming
conventions. The options available are documented in ``config.hpp``.

All configuration with config files is done for algorithm settings, which will be changed between
runs to find optimal parameters. All These options can be given in two ways: Specified in the 
config file or on command line. Command line options override config file options, allowing you to 
use a config file as a base but then change a few key options for between runs. 
Detailled explanations for all values can be found in the options.pdf file.

## Inputs

If not specified otherwise, values should be given in SI-units, e.g. meter for length, radians for
degree, seconds for time and Hertz for frequency. The output will be in SI-units as well.

Waveform data formatting
 - every line gives one time, amplitude pair
 - the order of values in a line is 1: time, 2: amplitude
 - Time unit should be seconds
 - Amplitude can be displacement, velocity or acceleration
 - values are separated by whitespace, multiple spaces are treated as one
 - the decimal separator is a dot, and scientific notation (1.2E-2) is understood
  

VelocityModel1D file formatting
 - lines prefixed with # are ignored -> comments
 - inline comments are not enabled
 - one line in the file describes one layer
 - whitespace is ignored
 - values are separated by commas, a dot is used as the decimal separator
 - no comma is required at the end of a line, the "\n" acts as a line separator
 - the dot can be omitted for integer values
 
 VelocityModel1D description
 - a VelocityModel is made up of one or multiple layers
 - a layer has multiple properties (p and s velocity, density)
 - the upper depth of the layer is inclusive, the lower depth is exclusive, this means that for a
 layer going from 5 to 10 km evaluating the layer at 5 km will yield the property of this layer,
 while evaluating it at 10 km will yield the property of the layer below
 - there are two type of layers: ConstantLayer and LinearLayer
   - a ConstantLayer has the same value for a property across its entire depth
   - a LinearLayers properties values changes with depth by a linear gradient given with a start
   value (top) and an end value (bottom)
 - a model can have either all linear or all constant layers
 - for a constant layer the order of values in a line is 
   depth_top, depth_bottom, p velocity, s velocity, density
 - for a linear layer the order of values in a line is 
   depth top, depth bottom, p velocity top, p velocity bottom, s velocity top, s velocity bottom, density top, density bottom

Project folder structure
 - Seismological data has to be stored in a certain structure. A project consists of data for a number of sources and receivers,
all in a common velocity model. 
 - The following structure is used for project folders:
   - Top level directory, specifies project name.
     - sources.txt, source file, specifies positions of sources in model.
     - receivers.txt, receiver file, specifies positions of receivers in model.
     - shotdata, directory containing seismograms.
       - source_0001, directory containing recordings from all seismograms for the shot.
         - receiver_001.txt, text file, format as given under Waveform data formatting, giving data recorded at this receiver.
         - The receivers index in the filename should be zero padded, and the receivers should have the same order as given in
         the receiver file, eg. index 1 in receiver file should correspond to receiver_001.txt. 
       - source_n
     - It is assumed that all seismograms share their timesteps, eg. all receiver files will have the same time values.
       It is also assumed, that this timestep is constant for one seismogram. 
     
   
## Literature/References

In the code, literature is often referred by a shorthand first author + year notation.
In case the title and full author list is not mentioned in the source code:  
Cerveny2001: Seismic ray theory - Cerveny 2001
Hill2001: Prestack Gaussian beam depth migration - Hill 1990
Fang2019: A fast and robust two-point ray tracing method in layered media with constant or linearly varying layer velocity - Fang, Chen 2019