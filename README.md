# ExternalFlowDataGeneration
Preliminary scripts generating external flow simulation examples for graph ML.

## Overview
Here we provide a repository for generation of large amounts of flow simulations with [fenics](https://fenicsproject.org/) and [oasis](https://github.com/mikaem/Oasis) solvers.
Currently cylinder geometries, ellipse geometries and channel nozzle geometries are supported. Geometries are generated and meshed with a [gmsh](https://gmsh.info/) kernel, and mesh formats are processed with [meshio](https://github.com/nschloe/meshio) for adaptation to FEniCS.

## Setup
To create a environment using conda:
```shell
conda update -n base -c defaults conda
conda env create --name fenicsproject --file=environment.yaml
conda activate fenicsproject
```

## Quick start

### Mesh generation
Mesh generation is completed in mesh.py. Change parameters in 
```python
if __name__ == "__main__":
```
part to create mesh of different sizes. Or write your own parameters and incorporate 
```python
from mesh import createFreeFlowMesh, createChannelFlowMesh
```
to use the [gmsh] kernel for your own mesh. Mesh file is stored by default in `io_operations`.

### Run simulations
For a free flow dataset, an oasis simulation can be executed by:
```shell
python NSfracStep.py problem=Cylinder meshname=ellipse_1 meshdir=io_operations T=10 dt=0.001 output_timeseries_as_vector=False folder=data
```
To run multiple simulations, refer to `simulation_series_oasis.py`. Simply run the script to perform simulation on all samples generated.
Run 'simulation_series.py' to perform the same simulation on a naive IPCS implementation by directly encoding all the functions into FEniCS.

### Construct dataset
