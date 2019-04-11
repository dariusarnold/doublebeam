# Double gaussian beam

TODO: short introduction

## Setup

TODO: add conda environment file

## Inputs

Waveform data formatting
 - every line gives one time, amplitude pair
 - the order of values in a line is time amplitude
 - TODO meaning of amplitude 
 - values are separated by whitespace, multiple spaces are ignored
 - the decimal separator is a dot, and scientific notation (1.2E-2) is understood
 - TODO unit of values
  

VelocityModel1D formatting
 - lines prefixed with # are ignored -> comments
 - in line comments are not enabled
 - whitespace is ignored
 - values are seperated by commas, a dot is used as the decimal seperator
 - the order of values is depth, p velocity, s velocity, density
 - the units are km, km/s, km/s, g/cm³
 - one line describes one layer where the depth gives the top of the layer
 - values are assumed to be constant between sample points
 - the depth marks the top of the layer, a layer extends downwards until another layer starts
 