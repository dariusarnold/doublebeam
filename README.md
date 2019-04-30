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
  

VelocityModel1D file formatting
 - lines prefixed with # are ignored -> comments
 - in line comments are not enabled
 - one line in the file describes one layer
 - whitespace is ignored
 - values are separated by commas, a dot is used as the decimal separator
 - no comma is required at the end of a line
 - the dot can be omitted for integer values
 
 VelocityModel1D description
 - a VelocityModel is made up of one or multiple layers
 - a layer has multiple properties (p and s velocity, density)
 - the upper depth of the layer is inclusive, the lower depth is exclusive
 - there are two type of layers: ConstantLayer and LinearLayer
   - a ConstantLayer has the same value for a property across its entire depth
   - a LinearLayers properties values changes with depth by a linear gradient given
     with a start value (top) and an end value (bottom)  
 - a model can have either all linear or all constant layers
 - for a constant layer the order of values in a line is 
   depth_top, depth_bottom, p velocity, s velocity, density
 - the units are km, km, km/s, km/s, g/cm³
 - for a linear layer the order of values in a line is 
   depth top, depth bottom, p velocity top, p velocity bottom, s velocity top, s velocity bottom, density top, density bottom 
 - the units for the top and bottom point are the same as in the constant layer case 
 