import os
from dolfin import *
from utils import *
import matplotlib.pyplot as plt

def plot_solution(mesh_name, mesh_directory, solution_directory, t):
    '''
    Plot fenics solution at different time steps.
    mesh_name: name of the mesh file series. e.g. 'circle_0'.
    mesh_directory: directory to mesh file.
    solution_directory: directory containing targeted solution file 'velocity_series' and 'pressure_series'.
    t: time step at which to plot the solution.
    '''
    mesh, mf_boundaries, association_table = import_mesh(prefix=mesh_name, subdomains=False, directory=mesh_directory)

    V = FunctionSpace(mesh, 'P', 2)
    Q = FunctionSpace(mesh, 'P', 1)
    u_  = Function(V)
    p_  = Function(Q)

    # velocity_series = TimeSeries(os.path.join(solution_directory, 'velocity_series'))
    # pressure_series = TimeSeries(os.path.join(solution_directory, 'pressure_series'))

    # velocity_series.retrieve(u_.vector(), t / 1000)
    # pressure_series.retrieve(p_.vector(), t / 1000)
    pressure_series = TimeSeries('oasis/pressure_series')
    pressure_series.retrieve(p_.vector(), t / 1000)

    velocity_series = TimeSeries('oasis/velocity_series')
    velocity_series.retrieve(u_.vector(), t / 1000)
    # velocity_series = XDMFFile('u_from_tstep_0.xdmf')
    # velocity_series.read(mesh)

    plot(u_)
    plt.show()

if __name__ == "__main__":
    simulation_type = 'ellipse'
    simulation_res = 'has'
    index = 3
    mesh_file = '{}_{}'.format(simulation_type, index)
    mesh_directory = "io_operations/data221129/{}/".format(simulation_res)
    solution_directory = 'navier_stokes_cylinder/sol221111_{}_{}'.format(index, simulation_res)
    t = 800

    plot_solution(mesh_file, mesh_directory, solution_directory, t)
