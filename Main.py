from classes.Grid import Grid
from classes.UniversalElement import UniversalElement
from classes.Element import Element
import classes.Functions_Const as consts
from classes.Node import Node
import numpy as np


def calculateElementField(element: Element, uniElem: UniversalElement) -> float:
    # claculate jacobian for each point in normal element
    jacobians = [uniElem.jacobian_matrix(element, n) for n in range(4)]
    # get weights for specified gauss points(here for 2. level there is only 1 as  )
    weights = []
    for y in range(2):
        for x in range(2):
            weights.append((consts.gauss_points[0][x][1], consts.gauss_points[0][y][1]))

    # calculating integral S = integral from -1 to 1 integral from -1 to 1 det[J] dKsi dEta using gauss quadrants
    sum = 0.

    for weight, jacobian in zip(weights, jacobians):
        sum += np.linalg.det(jacobian) * weight[0] * weight[1]

    return sum

if __name__ == "__main__":
    grid = Grid(4, 6, 0.5, 1.5)
    print(grid.getElement(6))
    grid.saveToJson("data.json")
    grid = Grid.LoadFromData("data.json")
    print(grid.getElement(6))
    uniElem = UniversalElement()
    print(uniElem.jacobian_matrix(grid.getElement(6), 0))
    print(UniversalElement.integral_gauss_3(lambda x: x**2, 1., 10.))
    #test for calculating element field
    nodes = [ Node(0, 1, 1, BC=True),
              Node(1, 3, 1, BC=True),
              Node(2, 3, 4, BC=True),
              Node(3, 1, 4, BC=True) ]
    element = Element(0, 4, nodes)
    print("Element field: {}".format(calculateElementField(element, uniElem)))


