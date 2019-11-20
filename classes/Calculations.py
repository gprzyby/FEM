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

def matrix_integral(jacobians: [], matrixes: []):
    # calculating integral of matrix
    to_ret = np.zeros(shape=(4, 4))
    weights = []
    for y in range(2):
        for x in range(2):
            weights.append((consts.gauss_points[0][x][1], consts.gauss_points[0][y][1]))

    for weight, jacobian, matrix in zip(weights, jacobians, matrixes):
        to_ret += matrix * weight[0] * weight[1] * np.linalg.det(jacobian)

    return to_ret


def distance_between_nodes(nodes: [Node]):
    nodes[0].x
    distance = (nodes[0].x - nodes[1].x) ** 2 + (nodes[0].y - nodes[1].y) ** 2
    return distance ** 0.5


# TODO: don't forget with integral with boundry conditions(not implemented yet)
def calculate_h_matrix(elem: Element, univElem: UniversalElement, conductivity: float, alpha: float):
    #creating jacobians for every integral point
    jacobians = [univElem.jacobian_matrix(elem, n) for n in range(4)]

    h = []
    for node_i in range(4):   #for 4 integral points
        # calculating dN_dx, dN_dy for every integral point
        dN_dx = []
        dN_dy = []
        for n_i in range(4):
            # making matrix of [[dN_dksi], [dN_deta]]
            temp = np.array([univElem.dN_dKsi[node_i][n_i], univElem.dN_dEta[node_i][n_i]]).reshape((2, 1))
            matrix_dn_dx_dy = np.dot(np.linalg.inv(jacobians[node_i]), temp)
            dN_dx.append(matrix_dn_dx_dy[0][0])
            dN_dy.append(matrix_dn_dx_dy[1][0])

        matrix_dN_dx = np.array([dN_dx])
        matrix_dN_dy = np.array([dN_dy])

        #TODO: check transpose
        temp = np.dot(np.transpose(matrix_dN_dx), matrix_dN_dx) + np.dot(np.transpose(matrix_dN_dy), matrix_dN_dy)
        h.append(temp)

    h_without_bc = matrix_integral(jacobians, h)

    # boundary conditions for newton law
    edges_with_bc = {}
    for i in range(4):
        if elem[i].BC and elem[(i + 1) % 4].BC:
            edges_with_bc[i] = [elem[i], elem[(i + 1) % 4]]

    # calculating jacobians for each edge (half of edge length)
    edges_jacobians = {}
    for key, nodes in edges_with_bc.items():
        edges_jacobians[key] = distance_between_nodes(nodes) / 2.

    # creating form functions for 1d
    N_dot_NT = np.zeros(shape=(2, 2))
    for integral_point in [consts.gauss_points[0][0][0], consts.gauss_points[0][1][0]]:
        N_1d = np.array([[consts.shape1d_fun[0](integral_point), consts.shape1d_fun[1](integral_point)]])
        temp = np.dot(np.transpose([N_1d]), N_1d)
        N_dot_NT += np.dot(np.transpose(N_1d), N_1d) * 1 * 1    # * detJ make later

    # combining matrices
    H_bc = np.zeros(shape=(4, 4))
    for bc_index in edges_with_bc:
        if bc_index != 3:       #in that case adding arrays is simple, just translation of N_dot_NT
            H_bc[bc_index:bc_index+2, bc_index:bc_index+2] += N_dot_NT * edges_jacobians[bc_index]
        else:
            H_bc[0:4:3, 0:4:3] += N_dot_NT * edges_jacobians
    H_bc *= alpha

    print(H_bc)

    return h_without_bc * conductivity + H_bc

def calculate_c_matrix(elem: Element, univ_elem: UniversalElement, density: float, specific_heat: float):
    # creating jacobians for every integral point
    jacobians = [univ_elem.jacobian_matrix(elem, n) for n in range(4)]
    n_matrix = univ_elem.N

    c_matrixes = []

    for point_n_matrix in n_matrix:
        tmp = np.dot(np.transpose(np.array([point_n_matrix])), np.array([point_n_matrix]))
        c_matrixes.append(tmp)

    to_ret = matrix_integral(jacobians, c_matrixes)

    return to_ret * density * specific_heat
