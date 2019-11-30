from classes.Grid import Grid
from classes.UniversalElement import UniversalElement
from classes.Element import Element
import classes.Functions_Const as consts
from classes.Node import Node
import numpy as np
import classes.Calculations as calc

if __name__ == "__main__":
    grid = Grid(4, 6, 0.5, 1.5)
    print(grid.getElement(6))
    grid.saveToJson("data.json")
    grid = Grid.LoadFromData("data.json")
    print(grid.getElement(6))
    uniElem = UniversalElement()
    print(uniElem.jacobian_matrix(grid.getElement(6), 0))
    print(calc.integral_gauss_3(lambda x: x**2, 1., 10.))
    #test for calculating element field
    nodes = [ Node(0, 0, 0, BC=True),
              Node(1, 0.025, 0, BC=True),
              Node(2, 0.025, 0.025, BC=True),
              Node(3, 0, 0.025, BC=False) ]
    element = Element(0, 4, nodes)
    print("Element field: {}".format(calc.calculateElementField(element, uniElem)))
    print(calc.calculate_h_matrix(element, uniElem, 30., 25))
    print(calc.calculate_c_matrix(element, uniElem, 7800, 700))
    print(calc.calculate_p_matrix(element, uniElem, 25, 1))
    test_grid = Grid(4, 4, 0.1, 0.1)
    print(calc.calculate_h_global(test_grid, uniElem, 25, 300))
    print(calc.calculate_c_global(test_grid, uniElem, 7800, 700))
    print(calc.calculate_p_global(test_grid, uniElem, 300, 1200))
    calc.simulate_heat_transfer(Grid(31, 31, 0.1, 0.1), 100, 100, 1, 1200, 300, 700, 25, 7800)

