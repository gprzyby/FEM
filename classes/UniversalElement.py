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

    @staticmethod
    def integral_gauss_3(function, x_beg, x_end):
        gauss_points = [gauss[0] for gauss in funConst.gauss_points[1]]
        gauss_weights = [gauss[1] for gauss in funConst.gauss_points[1]]
        det_jacobian = funConst.shape1d_dksi[0] * x_beg + funConst.shape1d_dksi[1] * x_end

        global_points = [funConst.shape1d_fun[0](point) * x_beg + funConst.shape1d_fun[1](point) * x_end for point in gauss_points]

        integral = sum([function(glob_point) * weight for glob_point, weight in zip(global_points, gauss_weights)]) * det_jacobian

        return integral