# glyco-analysis

### Introduction:

Glycogen_analysis.py is a python-blender API based script that performs analysis on a reconstructed module of glycogen data. The work presented here is part of the ongoing Neuro-Inspired High Performance Computing project, which is a collaboration between KAUST core labs and the BBP team at EPFL. This work is also considered part of a number of stages for implementing Neural Energetics from a visualization perspective.

#### Software and Data requirements:
To date, the addon has been maintaned and tested to work with blender 2.79b. For this, python3.5.3 is a requirement for the implementation of the DBSCAN algorithm and nearest neighbour search. The Python bundle should have site-packages support that includes:
- Numpy
- scipy
- sklearn
- matplotlib

The dataset should follow certain conventions such as names and objects-category-layers, this is mostly relates to Blender.
More details can be found in the [wiki][]
[wiki]:https://github.com/daniJb/glyco-analysis/wiki


Glycogen analysis has been tested on Linux centOS7, MacOSX Sierra and Windows10.

#### Addon features:
- hide/unhide group of objects based on selection.
- measurements as per granule / cluster.
- visualization of clusters with ellipsoids.
- calculate optimum clustering parameters for DBSCAN algorithm with interactive graph display.
- export measurements data to file.

