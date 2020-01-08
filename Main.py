from classes.Grid import Grid
from classes.UniversalElement import UniversalElement
from classes.Element import Element
import classes.Functions_Const as consts
from classes.Node import Node
import numpy as np
import classes.Calculations as calc
import sys
import json

"""
    uzyte warstwy scian:
    zwykla sciana nieocieplona: 15 mm tynku betonowego(ściślej mówiąć cementowo-wapienny)
                                25 cm muru ceglanego
    ściana ocieplona:   15 mm tynku betonowego(ściślej mówiąć cementowo-wapienny)
                        25 cm muru ceglanego  
                        5 cm styropianu
    ściana starego domu:    25 mm tynku betonowego
                            30 cm drewna
"""


def building_wall_simulation():
    print("Liczienie temperatur dla ściany nieocieplonej")
    layer_info = [[0.015, {"conductivity": 2, "density": 2500, "specific_heat": 1000}],
                  [0.26, {"conductivity": 0.5, "density": 1000, "specific_heat": 2100}]]
    boundary_info = {0: {"alpha": 0, "ambient_temp": 0},
                     1: {"alpha": 30, "ambient_temp": -10},
                     2: {"alpha": 0, "ambient_temp": 0},
                     3: {"alpha": 0, "ambient_temp": 0}}
    grid = Grid(31, 31, 0.265, 0.1)
    calc.simulate_heat_transfer(grid, 21, 172800, 30, layer_info, boundary_info)
    grid.saveToVTK("./grid_nieocieplonaSciana")

    print("Liczienie temperatur dla ściany starego domu")
    layer_info = [[0.035, {"conductivity": 2, "density": 2500, "specific_heat": 1000}],
                  [0.32, {"conductivity": 0.5, "density": 1000, "specific_heat": 2100}]]
    boundary_info = {0: {"alpha": 0, "ambient_temp": 0},
                     1: {"alpha": 30, "ambient_temp": -10},
                     2: {"alpha": 0, "ambient_temp": 0},
                     3: {"alpha": 0, "ambient_temp": 0}}
    grid = Grid(31, 31, 0.335, 0.1)
    calc.simulate_heat_transfer(grid, 21, 172800, 30, layer_info, boundary_info)
    grid.saveToVTK("./grid_staryDom")

    print("Liczienie temperatur dla ściany ocieplonej")
    layer_info = [[0.015, {"conductivity": 2, "density": 2500, "specific_heat": 1000}],
                  [0.26, {"conductivity": 0.5, "density": 1000, "specific_heat": 2100}],
                  [0.05, {"conductivity": 0.04, "density": 30, "specific_heat": 1500}]]
    boundary_info = {0: {"alpha": 0, "ambient_temp": 0},
                     1: {"alpha": 30, "ambient_temp": -10},
                     2: {"alpha": 0, "ambient_temp": 0},
                     3: {"alpha": 0, "ambient_temp": 0}}
    grid = Grid(31, 31, 0.325, 0.1)
    calc.simulate_heat_transfer(grid, 21, 172800, 30, layer_info, boundary_info)
    grid.saveToVTK("./grid_ocieplonaSciana")


def main():
    if (len(sys.argv) == 1):
        print("Run with parameters: python Main.py <simulation data in json path> <filename to save in paraview format>")
        exit(-1)
    uniElem = UniversalElement()
    file_path = sys.argv[1]
    try:
        grid = Grid.LoadFromData(file_path)
        with open(sys.argv[1]) as data_file:
            data_json = json.load(data_file)
            alpha = data_json["alpha"]
            ambient_temp = data_json["ambient_temp"]
            conductivity = data_json["conductivity"]
            specific_heat = data_json["specific_heat"]
            density = data_json["density"]
            simulation_time = data_json["simulation_time"]
            simulation_step = data_json["step_time"]
            initial_temp = data_json["initial_temp"]
            boundary_info = {0: {"alpha": alpha, "ambient_temp": ambient_temp},
                             1: {"alpha": alpha, "ambient_temp": ambient_temp},
                             2: {"alpha": alpha, "ambient_temp": ambient_temp},
                             3: {"alpha": alpha, "ambient_temp": ambient_temp}}
            # setting last layer bigger to eliminate float error
            layer_info = [[grid.get_spec().W + 1, {"conductivity": conductivity, "density": density, "specific_heat": specific_heat}]]
            calc.simulate_heat_transfer(grid, initial_temp, simulation_time, simulation_step, layer_info, boundary_info)
            #save to vtk
            if(len(sys.argv) == 3):
                grid.saveToVTK(sys.argv[2])

    except IOError:
        print("File not found!")


if __name__ == "__main__":
    main()


