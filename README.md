# Ray Casting on Triangle Meshes

## Overview

  Projection and animation of arbitrary number of objects given in obj format using a computer graphics technique called "Ray casting" to make them look realistic.
  
  Refer the documentation for more details https://ksskreddy.github.io/RayTracing-on-Triangle-Meshes/

 ## Specifications
  
  <b>Input</b> : Set of objects(in obj format), respective centers, scaling factors, colors, position of light and camera

  <b>Output</b>: The scene as seen from the camera in a 500×500 PPM image

 ### Other details

   * The objects are assumed to be in a 500×500×500 room centered at (250,250,−250)
   * So, the centers of the objects and their scaling factors must be chosen accordingly so that the objects don't go outside the room.
   * The specifications like room size, output image dimensions can be changed by changing the macro WALL_SIDE in the tracer.h file.
   * To turn off the effect of shadow of one object on another, the variable CrossShadows can be set to false.
   
## Results

### Scene with vase on table

  <p align="center"> <img align = "center" src= "https://github.com/ksskreddy/Graphics-Project/blob/master/assets/images/table.gif" title = "Scene with vase on table"></p>

### Scene with teapot and mug

  <p align="center"> <img align = "center" src= "https://github.com/ksskreddy/Graphics-Project/blob/master/assets/images/teapot.gif" title = "Scene with teapot and mug"></p>

### Scene with Spheres

  <p align="center"> <img align = "center" src= "https://github.com/ksskreddy/Graphics-Project/blob/master/assets/images/spheres.gif" title = "Scene with Spheres"></p>

### Scene with bunny

  <p align="center"> <img align = "center" src= "https://github.com/ksskreddy/Graphics-Project/blob/master/assets/images/bunny.gif" title = "Scene with bunny"></p>
 
### Scene with plane

  <p align="center"> <img align = "center" src= "https://github.com/ksskreddy/Graphics-Project/blob/master/assets/images/plane.gif" title = "Scene with plane"></p>


