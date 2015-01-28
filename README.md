# glyco-analysis

### Introduction:

Glycogen_analysis.py is a python-blender API based script that performs analysis on a reconstructed module of glycogen data. The work presented here is part of the ongoing Neuro-Inspired High Performance Computing project. It is a collaboration between KAUST core labs and the BBP team at EPFL. This work is also considered part of a number of stages for implementing Neural Energetics from a visualization perspective.

#### Software and Data requirements:
For the addon to run properly it requires python3.4.2 and blender2.72b. Python installation should have site-packages support that includes:
- Numpy
- scipy
- sklearn
- matplotlib

The dataset should follow certain conventions such as names and objects-category-layers, this is mostly relates to Blender.
More details can be found in the [wiki][]
[wiki]:https://github.com/daniJb/glyco-analysis/wiki

For better runtime performance, you need to set OMP_NUM_THREADS to 18. This is for an 8GB memory, 8 cores machine.
```
$ export OMP_NUM_THREADS=18
```

This was tested on Linux SL6, CentOS7 and MacOSX (Yosemite/Mountain Lion).

#### Addon features:
- hide/unhide group of objects based on selection.
- measurements as per granule / cluster.
- calculate optimum clustering parameters for DBSCAN algorithm with interactive graph display.
- save measurements data to file.

