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

def integral_gauss_3(function, x_beg, x_end):
    gauss_points = [gauss[0] for gauss in consts.gauss_points[1]]
    gauss_weights = [gauss[1] for gauss in consts.gauss_points[1]]
    det_jacobian = consts.shape1d_dksi[0] * x_beg + consts.shape1d_dksi[1] * x_end

    global_points = [consts.shape1d_fun[0](point) * x_beg + consts.shape1d_fun[1](point) * x_end for point in gauss_points]

    integral = sum([function(glob_point) * weight for glob_point, weight in zip(global_points, gauss_weights)]) * det_jacobian

    return integral

# TODO: don't forget with integral with boundry conditions(not implemented yet)
def calculate_h_matrix(elem: Element, univElem: UniversalElement):
    #creating jacobians for every integral point
    jacobians = [univElem.jacobian_matrix(elem, n) for n in range(4)]

    h = []
    for node_i in range(4):   #for 4 integral points
        # calculating dN_dx, dN_dy for every integral point
        dN_dx = []
        dN_dy = []
        for n_i in range(4):
            # making matrix of [[dN_dksi], [dN_deta]]
            temp = np.array([univElem.dN_dEta[node_i][n_i], univElem.dN_dKsi[node_i][n_i]]).reshape((2, 1))
            matrix_dn_dx_dy = np.dot(np.linalg.inv(jacobians[node_i]), temp)
            dN_dx.append(matrix_dn_dx_dy[0][0])
            dN_dy.append(matrix_dn_dx_dy[1][0])

        matrix_dN_dx = np.array([dN_dx])
        matrix_dN_dy = np.array([dN_dy])

        #TODO: check transpose
        temp = np.dot(np.transpose(matrix_dN_dx), matrix_dN_dx) + np.dot(np.transpose(matrix_dN_dy), matrix_dN_dy)
        h.append(temp)

    #calculating integral of matrix
    to_ret = np.zeros(shape=(4, 4))
    weights = []
    for y in range(2):
        for x in range(2):
            weights.append((consts.gauss_points[0][x][1], consts.gauss_points[0][y][1]))

    for weight, jacobian, h_in in zip(weights, jacobians, h):
        to_ret += h_in * weight[0] * weight[1] * np.linalg.det(jacobian)

    return to_ret

