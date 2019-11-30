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
    distance = (nodes[0].x - nodes[1].x) ** 2 + (nodes[0].y - nodes[1].y) ** 2
    return distance ** 0.5


def calculate_h_matrix(elem: Element, univElem: UniversalElement, conductivity: float, alpha: float):
    # creating jacobians for every integral point
    jacobians = [univElem.jacobian_matrix(elem, n) for n in range(4)]

    h = []
    for node_i in range(4):   # for 4 integral points
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

        temp = np.dot(np.transpose(matrix_dN_dx), matrix_dN_dx) + np.dot(np.transpose(matrix_dN_dy), matrix_dN_dy)
        h.append(temp)

    h_without_bc = matrix_integral(jacobians, h) * conductivity

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
        if bc_index != 3:       # in that case adding arrays is simple, just translation of N_dot_NT
            H_bc[bc_index:bc_index+2, bc_index:bc_index+2] += N_dot_NT * edges_jacobians[bc_index]
        else:
            H_bc[0:4:3, 0:4:3] += N_dot_NT * edges_jacobians[bc_index]
    H_bc *= alpha


    return h_without_bc + H_bc


def calculate_c_matrix(elem: Element, univ_elem: UniversalElement, density: float, specific_heat: float):
    # creating jacobians for every integral point
    jacobians = [univ_elem.jacobian_matrix(elem, n) for n in range(4)]
    n_matrix = univ_elem.N

    c_matrices = []

    for point_n_matrix in n_matrix:
        tmp = np.dot(np.transpose(np.array([point_n_matrix])), np.array([point_n_matrix]))
        c_matrices.append(tmp)

    to_ret = matrix_integral(jacobians, c_matrices)

    return to_ret * density * specific_heat


def calculate_p_matrix(elem: Element, univ_elem: UniversalElement, alpha: float, ambient_temp: float):
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
    N = np.zeros(shape=(2, 1))
    for integral_point in [consts.gauss_points[0][0][0], consts.gauss_points[0][1][0]]:
        N += np.array([[consts.shape1d_fun[0](integral_point), consts.shape1d_fun[1](integral_point)]]).reshape(2, 1)

    # combining matrices
    P = np.zeros(shape=(4, 1))
    for bc_index in edges_with_bc:
        if bc_index != 3:  # in that case adding arrays is simple, just translation of N_dot_NT
            P[bc_index:bc_index + 2, 0:1] += N * edges_jacobians[bc_index]
        else:
            P[0:4:3, 0:1] += N * edges_jacobians[bc_index]

    alpha *= -1
    P *= alpha * ambient_temp
    return P


def calculate_h_global(grid: Grid, univ_elem: UniversalElement, conductivity: float, alpha: float):
    global_h = np.zeros(shape=(grid.get_spec().mN, grid.get_spec().mN))

    for i in range(grid.get_spec().mE):
        element = grid.getElement(i + 1)
        local_h = calculate_h_matrix(element, univ_elem, conductivity, alpha)

        indexes = [element[i].id - 1 for i in range(4)]  # -1 to align with matrix indexing

        for local_row_index, global_row_index in zip(range(4), indexes):
            for local_column_index, global_column_index in zip(range(4), indexes):
                global_h[global_row_index, global_column_index] += local_h[local_row_index, local_column_index]

    return global_h


def calculate_c_global(grid: Grid, univ_elem: UniversalElement, density: float, specific_heat: float):
    global_c = np.zeros(shape=(grid.get_spec().mN, grid.get_spec().mN))

    for i in range(grid.get_spec().mE):
        element = grid.getElement(i + 1)
        local_c = calculate_c_matrix(element, univ_elem, density, specific_heat)

        indexes = [element[i].id - 1 for i in range(4)]  # -1 to align with matrix indexing

        for local_row_index, global_row_index in zip(range(4), indexes):
            for local_column_index, global_column_index in zip(range(4), indexes):
                global_c[global_row_index, global_column_index] += local_c[local_row_index, local_column_index]

    return global_c


def calculate_p_global(grid: Grid, univ_ele: UniversalElement, alpha: float, specific_heat: float):
    nodes_count = grid.get_spec().mN
    global_p = np.zeros(shape=(nodes_count, 1))

    for i in range(grid.get_spec().mE):
        element = grid.getElement(i + 1)
        local_p = calculate_p_matrix(element, univ_ele, alpha, specific_heat)

        indexes = [element[i].id - 1 for i in range(4)]  # -1 to align with matrix indexing

        for local_row_index, global_row_index in zip(range(4), indexes):
            global_p[global_row_index, 0] += local_p[local_row_index, 0]

    return global_p

def simulate_heat_transfer(grid: Grid, initial_temp: float, simulation_time: float, step_time: float, ambient_temp: float,
                           alpha: float, specific_heat: float, conductivity: float, density: float):
    nodes_temp = np.array([[initial_temp for _ in range(grid.get_spec().mN)]]).reshape((grid.get_spec().mN, 1))
    univ_elem = UniversalElement()
    steps = (int)(simulation_time/step_time)

    steps_times = [(step + 1) * step_time for step in range(steps)]

    for step in steps_times:
        global_h = calculate_h_global(grid, univ_elem, conductivity, alpha)
        global_c = calculate_c_global(grid, univ_elem, density, specific_heat)
        global_p = calculate_p_global(grid, univ_elem, alpha, ambient_temp)
        global_p *= -1
        global_c *= 1/step_time
        global_h += global_c
        global_p += np.dot(global_c, nodes_temp)
        nodes_temp = np.dot(np.linalg.inv(global_h), global_p)
        print("Step time: {}\tMin temp: {}\tMax temp: {}".format(step, np.amin(nodes_temp), np.amax(nodes_temp)))
