from classes.Node import Node
from classes.Element import Element
import json
import evtk.hl as vtk_export
import numpy as np

class Grid:

    # class used for serialization and deserialization to json
    class __grid_data:
        def __init__(self, H: float, W:float, mH: int, mW: int, mN: int, mE: int):
            self.H = H
            self.W = W
            self.mH = mH
            self.mW = mW
            self.mN = mN
            self.mE = mE

    def __init__(self, n_horizontal_nodes: int, n_vertical_nodes: int, gridWidth: float, gridHeight: float):
        # checking input
        if n_horizontal_nodes <= 1 or n_vertical_nodes <= 1:
            raise ValueError("Grid must have positive numbers or vertical and horizontal nodes")
        #initializing class fields
        nElements = (n_horizontal_nodes - 1) * (n_vertical_nodes - 1)
        nNodes = n_horizontal_nodes * n_vertical_nodes
        self.__data = self.__grid_data(gridHeight, gridWidth, n_vertical_nodes, n_horizontal_nodes, \
                                       nNodes, nElements)
        self.__nodes = []
        self.__elements = []
        element_width = gridWidth / (n_horizontal_nodes - 1)
        element_height = gridHeight / (n_vertical_nodes - 1)

        #creating list of nodes
        id = 1
        for x in range(n_horizontal_nodes):
            for y in range(n_vertical_nodes):
                #border control
                if (x == 0 or x == n_vertical_nodes - 1) or (y == 0 or y == n_horizontal_nodes - 1):
                    self.__nodes.append(Node(id, x * element_width, y * element_height, BC=True))
                else:
                    self.__nodes.append(Node(id, x * element_width, y * element_height, BC=False))
                id += 1

        #creating elements of 4 nodes
        elem_id = 1
        for x in range(n_horizontal_nodes - 1):
            for y in range(n_vertical_nodes - 1):
                id = y + x * n_vertical_nodes
                node_list = [self.__nodes[id], \
                             self.__nodes[n_vertical_nodes + id], \
                             self.__nodes[n_vertical_nodes + id + 1], \
                             self.__nodes[id + 1]]
                self.__elements.append(Element(elem_id, 4, node_list))
                elem_id += 1

    def getElement(self, id: int):
        return self.__elements[id - 1]

    def getNode(self, id: int):
        return self.__nodes[id - 1]

    def __str__(self):
        return "Grid"

    def get_spec(self):
        return self.__data

    @classmethod
    def LoadFromData(cls, path: str):
        with open(path, "r") as file:
            input = json.load(file)
            gridInfo = cls.__grid_data(**input)

        return cls(gridInfo.mW, gridInfo.mH, gridInfo.W, gridInfo.H)

    def saveToJson(self, path: str):
        with open(path, 'w') as file:
            json.dump(self.__data.__dict__, file)

    def saveToVTK(self, path: str):
        x_coords = [node.x for node in self.__nodes[::self.__data.mW]]
        y_coords = [node.y for node in self.__nodes[:self.__data.mH:]]
        temp = [node.value for node in self.__nodes]

        x = np.zeros((self.__data.mW, self.__data.mH, 1))
        y = np.zeros((self.__data.mW, self.__data.mH, 1))
        z = np.zeros((self.__data.mW, self.__data.mH, 1))

        for k in range(1):
            for j in range(self.__data.mH):
                for i in range(self.__data.mW):
                    x[i, j, k] = x_coords[i]
                    y[i, j, k] = y_coords[j]

        temp = np.asarray(temp).reshape((self.__data.mW, self.__data.mH, 1))

        vtk_export.gridToVTK(path, x, y, z, pointData={"Temperature": temp})