from classes.Node import Node

class Element(object):

    def __init__(self, id: int, nodes_num: int, node_list: [Node] = None):
        self.__nN = nodes_num
        self.__nodes = node_list
        self.__id = id

        if self.__nN <= 0:
            raise ValueError("Nodes number must be grater than zero")

        # checking if nodes number is same as len(__nodes)
        if len(self.__nodes) > self.__nN:
            raise ValueError("Parsed list has more elements than expected!")
        else:
            self.__nodes += [0.0 for _ in range(nodes_num - len(node_list))]

    @property
    def nodes_num(self):
        return self.__nN

    @property
    def id(self):
        return self.__id

    def __getitem__(self, item):
        return self.__nodes[item]

    #__setitem__ deleted for safety reasons(to have identity with grid class)

    def __str__(self):
        return "E id:{}, elements:".format(self.__id) + str(self.__nodes)

    def __repr__(self):
        return str(self)

    def points_coordinates_matrix(self):
        x = [node.x for node in self.__nodes]
        y = [node.y for node in self.__nodes]
        return x, y

