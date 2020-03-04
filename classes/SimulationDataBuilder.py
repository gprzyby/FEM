from enum import Enum


class SimulationDataBuilder:

    class Boundary(Enum):
        DOWN = 0
        RIGHT = 1
        UP = 2
        LEFT = 3

    def __init__(self):
        self.__layer_info = []
        default_boundary_data = {"alpha": 0., "ambient_temp": 0.}
        self.__boundary_cond = {i: default_boundary_data for i in range(4)}

    def add_layer(self, layer_thickness: float, **layer_data):
        self.__layer_info.append([layer_thickness, layer_data])
        return self

    def set_boundary(self, side: Boundary, **boundary_info):
        self.__boundary_cond[side.value] = boundary_info
        return self

    def build(self):
        layer_info = self.__layer_info
        boundary_info = self.__boundary_cond

        # clearing builder
        self.__init__()

        return layer_info, boundary_info

    def __repr__(self):
        return "SimulationDataBuilder(layer_data: {}, boundary_data: {}".format(self.__layer_info, self.__boundary_cond)
