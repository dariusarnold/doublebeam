# Double gaussian beam

TODO: short introduction

## Setup

TODO: add conda environment file

## Inputs

If not specified otherwise, values should be given in SI-units, e.g. meter for length, radians for
degree, seconds for time and kg for mass. The output will be in SI-units as well.

Waveform data formatting
 - every line gives one time, amplitude pair
 - TODO unit of time, amplitude
 - the order of values in a line is 1: time, 2: amplitude
 - TODO meaning of amplitude. Can displacement, velocity and acceleration be used?
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