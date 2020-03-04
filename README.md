# Finite element methods

Code for FEM projects in AGH

## Features

* calculating heat transfer for basic grid using FEM method 
* ability to change amount of integral points
* saving results to vtk format(compatible with paraview program)
* supports simulation for grid with different material layers

## Usage

To run program there must be json file with simulation data(thermal properties of material, simulation time etc.) Path to that file we put as first parameter of program. Second argument is optional and tells program where to save simualtion results(without extension)

Exmple run
```bash
python Main.py test_case1.json result_test1
``` 
	

