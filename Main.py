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
    nodes = [ Node(0, 1, 1, BC=True),
              Node(1, 3, 1, BC=True),
              Node(2, 3, 4, BC=True),
              Node(3, 1, 4, BC=True) ]
    element = Element(0, 4, nodes)
    print("Element field: {}".format(calc.calculateElementField(element, uniElem)))
    print(calc.calculate_h_matrix(element, uniElem))



