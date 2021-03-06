import numpy as np
import classes.Functions_Const as funConst
from classes.Element import Element

class UniversalElement:

    def __init__(self, gauss_level: int = 2):
        """
        Creates class with defined gauss level.
        :param gauss_level: number larger than zero
        """

        #for 2d element
        self.__gauss_level = gauss_level
        self.__dN_dksi_matrix = np.zeros((gauss_level ** 2, 4))    # for 4 nodes in element
        self.__dN_deta_matrix = np.zeros((gauss_level ** 2, 4))
        self.__N_matrix = np.zeros((gauss_level ** 2, 4))
        self.__weights = []
        self.__create_matrices()

    def __create_matrices(self):
        #crating integral points and weights
        integral_points = []
        for ksi_ip, ksi_weight in funConst.gauss_points[self.__gauss_level - 1]:
            for eta_ip, eta_weight in funConst.gauss_points[self.__gauss_level - 1]:
                integral_points.append((ksi_ip, eta_ip))
                self.__weights.append((ksi_weight, eta_weight))

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

    def jacobian_matrix(self, element: Element, integral_point_id):
        jacobian = np.zeros(shape=(2, 2))   # creating matrix 2x2 from NumPy
        x_arr, y_arr = element.points_coordinates_vector()  # creating list of nodes coordinates
        jacobian.itemset((0, 0), np.dot(self.__dN_dksi_matrix[integral_point_id], np.array([x_arr]).transpose())[0])    # dx_dksi
        jacobian.itemset((1, 0), np.dot(self.__dN_deta_matrix[integral_point_id], np.array([x_arr]).transpose())[0])    # dx_deta
        jacobian.itemset((1, 1), np.dot(self.__dN_deta_matrix[integral_point_id], np.array([y_arr]).transpose())[0])    # dy_deta
        jacobian.itemset((0, 1), np.dot(self.__dN_dksi_matrix[integral_point_id], np.array([y_arr]).transpose())[0])    # dy_dksi
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

    @property
    def Weights(self):
        return self.__weights

    @property
    def IP_count(self):
        return self.__gauss_level ** 2

    @property
    def Gauss_level(self):
        return self.__gauss_level

