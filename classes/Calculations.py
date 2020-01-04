from classes.Grid import Grid
from classes.UniversalElement import UniversalElement
from classes.Element import Element
import classes.Functions_Const as consts
from classes.Node import Node
import numpy as np


def __calculate_property_for_integral_points(elem: Element, univ_el: UniversalElement, layer_info: dict,
                                             property_name: str) -> list:
    property_for_ip = []
    pos_x_for_ip = []
    x_pos, _ = elem.points_coordinates_vector()

    # calculating global pos of integral points(here only x coord)
    for shape_fun_vector in univ_el.N:
        pos_x_for_ip.append(np.dot(shape_fun_vector, np.asarray([x_pos]).reshape(len(univ_el.N), 1))[0])

    for integral_id in range(len(univ_el.N)):
        act_layer_border = 0.
        for i, [layer_thickness, layer_properties] in enumerate(layer_info):
            if layer_thickness + act_layer_border >= pos_x_for_ip[integral_id]:  # pos of first ip (most to left)
                property_for_ip.append(layer_properties[property_name])
                break
            else:  # element is in other layer, go next
                act_layer_border += layer_thickness

    return property_for_ip


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

    global_points = [consts.shape1d_fun[0](point) * x_beg + consts.shape1d_fun[1](point) * x_end for point in
                     gauss_points]

    integral = sum(
        [function(glob_point) * weight for glob_point, weight in zip(global_points, gauss_weights)]) * det_jacobian

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


def calculate_h_boundary_matrix(elem: Element, boundaries_info: dict):
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
        N_dot_NT += np.dot(np.transpose(N_1d), N_1d) * 1 * 1  # * detJ make later

    # combining matrices
    H_bc = np.zeros(shape=(4, 4))
    for bc_index in edges_with_bc:
        if bc_index != 3:  # in that case adding arrays is simple, just translation of N_dot_NT
            H_bc[bc_index:bc_index + 2, bc_index:bc_index + 2] += N_dot_NT * edges_jacobians[bc_index] * \
                                                                  boundaries_info[bc_index]["alpha"]
        else:
            H_bc[0:4:3, 0:4:3] += N_dot_NT * edges_jacobians[bc_index] * boundaries_info[bc_index]["alpha"]

    return H_bc


def calculate_h_matrix(elem: Element, univElem: UniversalElement, layer_info: dict, boundaries_info: dict):
    # creating jacobians for every integral point
    jacobians = [univElem.jacobian_matrix(elem, n) for n in range(4)]

    # calculating conductivity and alpha for each integral point
    conductivity_for_ip = __calculate_property_for_integral_points(elem, univElem, layer_info, "conductivity")

    h = []
    for ip_number, node_i in enumerate(range(4)):  # for 4 integral points
        # calculating dN_dx, dN_dy for every integral point
        dN_dx = []
        dN_dy = []
        for n_i in range(4):  # for 4 shape functions
            # making matrix of [[dN_dksi], [dN_deta]]
            temp = np.array([univElem.dN_dKsi[node_i][n_i], univElem.dN_dEta[node_i][n_i]]).reshape((2, 1))
            matrix_dn_dx_dy = np.dot(np.linalg.inv(jacobians[node_i]), temp)
            dN_dx.append(matrix_dn_dx_dy[0][0])
            dN_dy.append(matrix_dn_dx_dy[1][0])

        matrix_dN_dx = np.array([dN_dx])
        matrix_dN_dy = np.array([dN_dy])

        temp = np.dot(np.transpose(matrix_dN_dx), matrix_dN_dx) + np.dot(np.transpose(matrix_dN_dy), matrix_dN_dy)
        temp *= conductivity_for_ip[ip_number]
        h.append(temp)

    h_without_bc = matrix_integral(jacobians, h)

    return h_without_bc + calculate_h_boundary_matrix(elem, boundaries_info)


def calculate_c_matrix(elem: Element, univ_elem: UniversalElement, layer_info: dict):
    # creating jacobians for every integral point
    jacobians = [univ_elem.jacobian_matrix(elem, n) for n in range(4)]
    n_matrix = univ_elem.N
    density_for_ip = __calculate_property_for_integral_points(elem, univ_elem, layer_info, "density")
    specific_heat_for_ip = __calculate_property_for_integral_points(elem, univ_elem, layer_info, "specific_heat")

    c_matrices = []

    for density_in_ip, specific_heat_in_ip, point_n_matrix in zip(density_for_ip, specific_heat_for_ip, n_matrix):
        tmp = np.dot(np.transpose(np.array([point_n_matrix])), np.array([point_n_matrix]))
        tmp *= density_in_ip * specific_heat_in_ip
        c_matrices.append(tmp)

    to_ret = matrix_integral(jacobians, c_matrices)

    return to_ret


def calculate_p_matrix(elem: Element, univ_elem: UniversalElement, boundaries_info: dict):
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
            P[bc_index:bc_index + 2, 0:1] += N * edges_jacobians[bc_index] * boundaries_info[bc_index]["alpha"] * \
                                             boundaries_info[bc_index]["ambient_temp"]
        else:
            P[0:4:3, 0:1] += N * edges_jacobians[bc_index] * boundaries_info[bc_index]["alpha"] * \
                             boundaries_info[bc_index]["ambient_temp"]

    P *= -1
    return P


def calculate_h_global(grid: Grid, univ_elem: UniversalElement, layer_info: dict, boundaries_info: dict):
    global_h = np.zeros(shape=(grid.get_spec().mN, grid.get_spec().mN))

    for i in range(grid.get_spec().mE):
        element = grid.getElement(i + 1)
        local_h = calculate_h_matrix(element, univ_elem, layer_info, boundaries_info)

        indexes = [element[i].id - 1 for i in range(4)]  # -1 to align with matrix indexing

        for local_row_index, global_row_index in zip(range(4), indexes):
            for local_column_index, global_column_index in zip(range(4), indexes):
                global_h[global_row_index, global_column_index] += local_h[local_row_index, local_column_index]

    return global_h


def calculate_c_global(grid: Grid, univ_elem: UniversalElement, layer_info: dict):
    global_c = np.zeros(shape=(grid.get_spec().mN, grid.get_spec().mN))

    for i in range(grid.get_spec().mE):
        element = grid.getElement(i + 1)
        local_c = calculate_c_matrix(element, univ_elem, layer_info)

        indexes = [element[i].id - 1 for i in range(4)]  # -1 to align with matrix indexing

        for local_row_index, global_row_index in zip(range(4), indexes):
            for local_column_index, global_column_index in zip(range(4), indexes):
                global_c[global_row_index, global_column_index] += local_c[local_row_index, local_column_index]

    return global_c


def calculate_p_global(grid: Grid, univ_ele: UniversalElement, boundaries_info: dict):
    nodes_count = grid.get_spec().mN
    global_p = np.zeros(shape=(nodes_count, 1))

    for i in range(grid.get_spec().mE):
        element = grid.getElement(i + 1)
        local_p = calculate_p_matrix(element, univ_ele, boundaries_info)

        indexes = [element[i].id - 1 for i in range(4)]  # -1 to align with matrix indexing

        for local_row_index, global_row_index in zip(range(4), indexes):
            global_p[global_row_index, 0] += local_p[local_row_index, 0]

    return global_p


def simulate_heat_transfer(grid: Grid, initial_temp: float, simulation_time: float, step_time: float, layer_info: dict,
                           boundaries_info: dict):
    nodes_temp = np.array([[initial_temp for _ in range(grid.get_spec().mN)]]).reshape((grid.get_spec().mN, 1))
    univ_elem = UniversalElement()
    steps = (int)(simulation_time / step_time)

    steps_times = [(step + 1) * step_time for step in range(steps)]

    for step in steps_times:
        global_h = calculate_h_global(grid, univ_elem, layer_info, boundaries_info)
        global_c = calculate_c_global(grid, univ_elem, layer_info)
        global_p = calculate_p_global(grid, univ_elem, boundaries_info)
        global_p *= -1
        global_c *= 1 / step_time
        global_h += global_c
        global_p += np.dot(global_c, nodes_temp)
        nodes_temp = np.dot(np.linalg.inv(global_h), global_p)
        print("Step time: {}\tMin temp: {}\tMax temp: {}".format(step, np.amin(nodes_temp), np.amax(nodes_temp)))

    # upgrade grid with temperatures
    for node_id, temp in enumerate(nodes_temp):
        grid.getNode(node_id + 1).value = temp[0]

    return nodes_temp
