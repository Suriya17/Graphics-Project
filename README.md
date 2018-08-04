# Ray Casting on Triangle Meshes

## Overview

  Projection and animation of arbitrary number of objects given in obj format using a computer graphics technique called "Ray casting" to make them look realistic.

 ## Specifications
  
  <b>Input</b> : Set of objects(in obj format), respective centers, scaling factors, colors, position of light and camera

  <b>Output</b>: The scene as seen from the camera in a 500×500 PPM image

 ### Other details

   The objects are assumed to be in a 500×500×500 room centered at (250,250,−250)
   So, the centers of the objects and their scaling factors must be chosen accordingly so that the objects don't go outside the room.
   The specifications like room size, output image dimensions can be changed by changing the macro WALL_SIDE in the tracer.h file.
To turn off the effect of shadow of one object on another, the variable CrossShadows can be set to false.
