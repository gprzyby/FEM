import numpy as np
import classes.Functions_Const as funConst
from classes.Element import Element

# TODO: try to make universal element for specified gaussian level(2., 3.)
class UniversalElement:

    def __init__(self):
        #for 2d element
        self.__dN_dksi_matrix = np.zeros((4, 4))    # for 4 nodes in element
        self.__dN_deta_matrix = np.zeros((4, 4))
        self.__N_matrix = np.zeros((4, 4))
        self.__create_matrices()

    def __create_matrices(self):
        integral_points = ((funConst.gauss_points[0][0][0], funConst.gauss_points[0][0][0]), \
                        (funConst.gauss_points[0][1][0], funConst.gauss_points[0][0][0]), \
                        (funConst.gauss_points[0][1][0], funConst.gauss_points[0][1][0]), \
                        (funConst.gauss_points[0][0][0], funConst.gauss_points[0][1][0]))

        #creating dN_dksi_matrix
        for point, array in enumerate(self.__dN_dksi_matrix):
            eta: float = integral_points[point][1]
            ksi: float = integral_points[point][0]
            for shape_fun_num, element in enumerate(array):
                np.put(array, shape_fun_num, funConst.shape2d_dksi[shape_fun_num](ksi=ksi, eta=eta))

        #creating dN_deta_matrix
        for point, array in enumerate(self.__dN_deta_matrix):
            eta: float = integral_points[point][1]
            ksi: float = integral_points[point][0]
            for shape_fun_num, element in enumerate(array):
                np.put(array, shape_fun_num, funConst.shape2d_deta[shape_fun_num](ksi=ksi, eta=eta))

        #creating N_matrix
        for point, array in enumerate(self.__N_matrix):
            eta: float = integral_points[point][1]
            ksi: float = integral_points[point][0]
            for shape_fun_num, element in enumerate(array):
                np.put(array, shape_fun_num, funConst.shape2d_fun[shape_fun_num](ksi=ksi, eta=eta))

    def jacobian_matrix(self, element: Element, element_num):
        #creating matrices from element points
        jacobian = np.zeros(shape=(2, 2))
        x_arr, y_arr = element.points_coordinates_matrix()
        print(np.dot(self.__dN_dksi_matrix[0],np.array([x_arr]).transpose())[0])
        jacobian.itemset((0, 0), np.dot(self.__dN_dksi_matrix[element_num],np.array([x_arr]).transpose())[0])    # dx_dksi
        jacobian.itemset((1, 0), np.dot(self.__dN_deta_matrix[element_num],np.array([x_arr]).transpose())[0])    # dx_deta
        jacobian.itemset((1, 1), np.dot(self.__dN_deta_matrix[element_num],np.array([y_arr]).transpose())[0])    # dy_deta
        jacobian.itemset((0, 1), np.dot(self.__dN_dksi_matrix[element_num],np.array([y_arr]).transpose())[0])    # dy_dksi
        return jacobian

    @property
    def dN_dKsi(self):
        return self.__dN_dksi_matrix

    @property
    def dN_dEta(self):
        return self.__dN_deta_matrix

    @property
    def N(self):
        return self.__N_matrix

