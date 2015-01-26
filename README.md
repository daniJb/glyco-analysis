# glyco-analysis

### Introduction:

Glycogen_analysis.py is a python-blender API based script that performs analysis on a reconstructed module of glycogen data. The work presented here is part of the ongoing Neuro-Inspired High Performance Computing project. It is a collaboration between KAUST core labs and the BBP team at EPFL. This work is also considered part of a number of stages for implementing Neural Energetics from a visualization perspective.

#### Requirements:
For the addon to run properly it requires python3.4.2 and blender2.72b. Python installation should have site-packages support that includes:
- Numpy
- scipy
- sklearn
- matplotlib

More installation details can be found here [link!]

For better runtime performance, you need to set OMP_NUM_THREADS to 18.
```
$ export OMP_NUM_THREADS=18
```

This was tested on Linux SL6, CentOS7 and MacOSX (Yusamite/Mountain Lion).

#### Addon features:
- hide/unhide group of objects based on selection.
- measurements as per granule / cluster.
- calculate optimum clustering parameters for DBSCAN algorithm with interactive graph display.
- save measurements data to file.

