# 🧬 Bio-Printing Trajectory Generation (SLERP & LERP)

This repository contains a spatial trajectory planning algorithm developed for bio-printing applications. The code processes 3D models in STL format, computes an optimized path over the surface, and ensures smooth orientation transitions of the printhead using Linear Interpolation (LERP) and Spherical Linear Interpolation (SLERP).

## Project Overview

In bio-printing over complex 3D geometries, it is essential to control not only the spatial position of the printhead \((X, Y, Z)\), but also its orientation relative to the surface. This project addresses that problem by generating a continuous trajectory that keeps the tool aligned with the local surface normals of the mesh.

The main workflow is:

1. **Load the STL mesh** and extract all triangular faces.
2. **Compute triangle centers and normal vectors** for the mesh.
3. **Convert normals into quaternions** so the tool orientation can be represented robustly in 3D space.
4. **Sort the points using a KDTree-based nearest-neighbor strategy** to build a more continuous path.
5. **Interpolate positions with LERP** and orientations with SLERP** to obtain a smooth, denser trajectory for execution or visualization.

## Requirements
This project uses Python and the following libraries:
```bash
pip install numpy matplotlib numpy-stl scipy
```
## Usage

1. Place your STL file in the root directory of the project.
2. Make sure the script points to the correct STL filename.
3. Run the program:
```bash
python slerp.py
```

## Visual Results
1. Original Mesh
![alt text](image.png)