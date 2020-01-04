from classes.Grid import Grid
from classes.UniversalElement import UniversalElement
from classes.Element import Element
import classes.Functions_Const as consts
from classes.Node import Node
import numpy as np
import classes.Calculations as calc

if __name__ == "__main__":
    uniElem = UniversalElement()
    alpha = 300
    ambient_temp = 1200
    layer_info = [[0.2, {"conductivity": 35, "density": 7800, "specific_heat": 700}]]
    boundary_info = {0: {"alpha": 0, "ambient_temp": 0},
                     1: {"alpha": 300, "ambient_temp": 1200},
                     2: {"alpha": 0, "ambient_temp": 0},
                     3: {"alpha": 300, "ambient_temp": 1200}}
    grid = Grid(4, 4, 0.1, 0.1)
    calc.simulate_heat_transfer(grid, 100, 500, 50, layer_info, boundary_info)
    grid.saveToVTK("./grid_paraviewNew")

