from classes.Grid import Grid
from classes.UniversalElement import UniversalElement
from classes.SimulationDataBuilder import SimulationDataBuilder
import classes.Calculations as calc
import sys
import json


def building_wall_simulation():
    univ_elem = UniversalElement()
    simulation_data_builder = SimulationDataBuilder()

    print("Liczienie temperatur dla ściany nieocieplonej")

    simulation_data_builder.add_layer(0.015, conductivity=2, density=2500, specific_heat=1000) \
        .add_layer(0.26, conductivity=0.5, density=1000, specific_heat=2100) \
        .set_boundary(SimulationDataBuilder.Boundary.RIGHT,
                      alpha=30, ambient_temp=-10)

    layer_info, boundary_info = simulation_data_builder.build()
    grid = Grid(31, 31, 0.265, 0.1)
    calc.simulate_heat_transfer(grid, univ_elem, 21, 17280, 30, layer_info, boundary_info)
    grid.saveToVTK("./grid_nieocieplonaSciana")

    print("Liczienie temperatur dla ściany starego domu")
    simulation_data_builder.add_layer(0.035, conductivity=2, density=2500, specific_heat=1000) \
        .add_layer(0.32, conductivity=0.5, density=1000, specific_heat=2100) \
        .set_boundary(SimulationDataBuilder.Boundary.RIGHT,
                      alpha=30, ambient_temp=-10)

    layer_info, boundary_info = simulation_data_builder.build()
    grid = Grid(31, 31, 0.335, 0.1)
    calc.simulate_heat_transfer(grid, univ_elem, 21, 17280, 30, layer_info, boundary_info)
    grid.saveToVTK("./grid_staryDom")

    print("Liczienie temperatur dla ściany ocieplonej")

    simulation_data_builder.add_layer(0.015, conductivity=2, density=2500, specific_heat=1000) \
        .add_layer(0.26, conductivity=0.5, density=1000, specific_heat=2100) \
        .add_layer(0.05, conductivity=0.04, density=30, specific_heat=1500) \
        .set_boundary(SimulationDataBuilder.Boundary.RIGHT,
                      alpha=30, ambient_temp=-10)

    layer_info, boundary_info = simulation_data_builder.build()
    grid = Grid(31, 31, 0.325, 0.1)
    calc.simulate_heat_transfer(grid, univ_elem, 21, 17280, 30, layer_info, boundary_info)
    grid.saveToVTK("./grid_ocieplonaSciana")


def main():
    if (len(sys.argv) == 1):
        print(
            "Run with parameters: python Main.py <simulation data in json path> <filename to save in paraview format>")
        exit(-1)
    file_path = sys.argv[1]
    try:
        grid = Grid.LoadFromData(file_path)
        with open(sys.argv[1]) as data_file:
            data_json = json.load(data_file)
            univ_elem = UniversalElement(4)
            calc.isotropic_heat_transfer_simulation(grid, univ_elem, **data_json)

            # save to vtk
            if (len(sys.argv) == 3):
                grid.saveToVTK(sys.argv[2])

    except IOError:
        print("File not found!")


if __name__ == "__main__":
    main()
