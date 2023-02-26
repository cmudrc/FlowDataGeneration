import os
from dolfin import *
from fenics import *
from utils import *
import matplotlib.pyplot as plt
from utils import *


def mesh2string(mesh):
    import matplotlib.tri as tri
    xy = mesh.coordinates()
    return tri.Triangulation(xy[:, 0], xy[:, 1], mesh.cells())


def plot_solution(mesh_name, mesh_directory, solution_directory, t):
    '''
    Plot fenics solution at different time steps.
    Input:
        mesh_name: name of the mesh file series. e.g. 'circle_0'.
        mesh_directory: directory to mesh file.
        solution_directory: directory containing targeted solution file 'velocity_series' and 'pressure_series'.
        t: time step at which to plot the solution.
    Ouput:
        None
    '''
    mesh, mf_boundaries, association_table = import_mesh(prefix=mesh_name, subdomains=False, directory=mesh_directory)

    V = FunctionSpace(mesh, 'P', 2)
    Q = FunctionSpace(mesh, 'P', 1)
    u0_  = Function(V)
    u1_  = Function(V)
    p_  = Function(Q)

    velocity_series_0 = TimeSeries(os.path.join(solution_directory, 'u0_from_tstep_0'))
    velocity_series_1 = TimeSeries(os.path.join(solution_directory, 'u1_from_tstep_0'))
    pressure_series = TimeSeries(os.path.join(solution_directory, 'p_from_tstep_0'))

    velocity_series_0.retrieve(u0_.vector(), t / 1000)
    velocity_series_1.retrieve(u0_.vector(), t / 1000)
    # pressure_series.retrieve(p_.vector(), t / 1000)
    # pressure_series = TimeSeries('oasis/pressure_series')
    pressure_series.retrieve(p_.vector(), t / 1000)
    X = mesh.coordinates()
    X = [X[:, i] for i in range(2)]
    possible_idx = []
    for i in range(len(X[0])):
        if X[0][i] > 9 and X[0][i] < 10 and X[1][i] > 5 and X[1][i] < 7:
            possible_idx.append(i)

    w0 = u0_.compute_vertex_values(mesh)
    # w1 = u1_.compute_vertex_values(mesh)
    w1 = u1_.compute_vertex_values(mesh)
    p1 = p_.compute_vertex_values(mesh)

    C_has = w0 ** 2
    C_has = np.sqrt(C_has)
    U = np.array([w0])
    args = X + U + [C_has]
    # plt.quiver(*args, scale=0.1)
    plot(u0_, mode='glyphs', scale=0.1, title='velocity')
    # velocity_series = TimeSeries('oasis/velocity_series')
    # velocity_series_0.retrieve(u_.vector(), t / 1000)
    # velocity_series = XDMFFile('u_from_tstep_0.xdmf')
    # velocity_series.read(mesh)
    # combine u0 and u1

    plot(u0_, mode='glyphs', scale=0.1, title='velocity')
    plt.show()


def plot_prediction(mesh_name, mesh_directiory, solution_directory):
    """
    Plot prediction from FlowML models
    Input:
        mesh_name: name of the mesh file series. e.g. 'circle_0'.
        mesh_directory: directory to mesh file.
        solution_directory: directory of numpy array containing the prediction, las and has.
        t: time step at which to plot the solution.
    """

    mesh, mf_boundaries, association_table = import_mesh(prefix=mesh_name, subdomains=False, directory=mesh_directory)

    gdim = mesh.geometry().dim()
    nv = mesh.num_vertices()

    prediction = np.load(solution_directory)
    u_has_star = prediction['u_has_star']
    p_has_star = prediction['p_has_star']
    u_has = prediction['u_has']
    p_has = prediction['p_has']
    u_las = prediction['u_las']
    p_las = prediction['p_las']

    X = mesh.coordinates()
    X = [X[:, i] for i in range(gdim)]

    C_has_star = u_has_star[0] ** 2
    for i in range(1, gdim):
        C_has_star += u_has_star[i] ** 2
    C_has_star = np.sqrt(C_has_star)

    # args = X + u_has_star + [C]
    # plt.figure()
    # ax.quiver(*args)
    # plt.title('velocity prediction at circle_109_909_model_325')

    C_has = u_has[0] ** 2
    for i in range(1, gdim):
        C_has += u_has[i] ** 2
    C_has = np.sqrt(C_has)

    # args = X + u_has + [C]
    # plt.figure()
    # ax.quiver(*args)
    # plt.title('velocity ground truth at circle_109_909_model_325')

    C_las = u_las[0] ** 2
    for i in range(1, gdim):
        C_las += u_las[i] ** 2
    C_las = np.sqrt(C_las)

    # calculate mean squared error between las, has and has_star, has
    mse_star = np.mean((u_has_star - u_has) ** 2)
    mse_star_ux = np.mean((u_has_star[0] - u_has[0]) ** 2)
    mse_star_uy = np.mean((u_has_star[1] - u_has[1]) ** 2)
    mse_star_p = np.mean((p_has_star - p_has) ** 2)
    print('mse between has and has_star ux: {}'.format(mse_star_ux))
    print('mse between has and has_star uy: {}'.format(mse_star_uy))
    print('mse between has and has_star p: {}'.format(mse_star_p))
    mse_las = np.mean((u_las - u_has) ** 2)
    mse_las_ux = np.mean((u_las[0] - u_has[0]) ** 2)
    mse_las_uy = np.mean((u_las[1] - u_has[1]) ** 2)
    mse_las_p = np.mean((p_las - p_has) ** 2)
    print('mse between las and has ux: {}'.format(mse_las_ux))
    print('mse between las and has uy: {}'.format(mse_las_uy))
    print('mse between las and has p: {}'.format(mse_las_p))
    relative_mse = (mse_star - mse_las) / mse_las
    print('mse between las and has: {}'.format(mse_las))
    print('mse between has and has_star: {}'.format(mse_star))

    # calculate maximum error between las, has and has_star, has
    max_error_star = np.max(np.abs(u_has_star - u_has))
    max_error_star_ux = np.max(np.abs(u_has_star[0] - u_has[0]))
    max_error_star_uy = np.max(np.abs(u_has_star[1] - u_has[1]))
    max_error_star_p = np.max(np.abs(p_has_star - p_has))
    print('max error between has and has_star ux: {}'.format(max_error_star_ux))
    print('max error between has and has_star uy: {}'.format(max_error_star_uy))
    print('max error between has and has_star p: {}'.format(max_error_star_p))
    max_error_las = np.max(np.abs(u_las - u_has))
    max_error_las_ux = np.max(np.abs(u_las[0] - u_has[0]))
    max_error_las_uy = np.max(np.abs(u_las[1] - u_has[1]))
    max_error_las_p = np.max(np.abs(p_las - p_has))
    print('max error between las and has ux: {}'.format(max_error_las_ux))
    print('max error between las and has uy: {}'.format(max_error_las_uy))
    print('max error between las and has p: {}'.format(max_error_las_p))
    relative_max_error = (max_error_star - max_error_las) / max_error_las
    print('max error between las and has: {}'.format(max_error_las))
    print('max error between has and has_star: {}'.format(max_error_star))

    # args = X + u_las + [C]
    # plt.figure()
    # ax.quiver(*args)
    # plt.title('velocity las at circle_109_909_model_325')

    # plt.figure()
    # plt.tricontourf(mesh2string(mesh), p_has_star, 40)
    # plt.axes('off')
    # plt.colorbar(label='pressure (Pa)')

    # plt.figure()
    # plt.tricontourf(mesh2string(mesh), p_has, 40)
    # plt.axes('off')
    # plt.colorbar(label='pressure (Pa)')

    # plt.figure()
    # plt.tricontourf(mesh2string(mesh), p_las, 40)
    # plt.axes('off')
    # plt.colorbar(label='pressure (Pa)')

    # plt.figure()
    # plt.tricontourf(mesh2string(mesh), p_has_star - p_has, 40)
    # plt.axes('off')
    # plt.colorbar(label='pressure (Pa)')

    # plt.figure()
    # plt.tricontourf(mesh2string(mesh), C_has_star, 40)
    # plt.axes('off')
    # plt.colorbar(label='velocity (m/s)')

    # plt.figure()
    # plt.tricontourf(mesh2string(mesh), C_has, 40)
    # plt.axes('off')
    # plt.colorbar(label='velocity (m/s)')

    # plt.figure()
    # plt.tricontourf(mesh2string(mesh), C_las, 40)
    # plt.axes('off')
    # plt.colorbar(label='velocity (m/s)')

    # plt.figure()
    # plt.tricontourf(mesh2string(mesh), C_has_star - C_has, 40)
    # plt.axes('off')
    # plt.colorbar(label='velocity (m/s)')
    # plt.show()


def plot_mesh(mesh_name, mesh_directory):
    """
    Plot mesh
    Input:
        mesh_name: name of the mesh file series. e.g. 'circle_0'.
        mesh_directory: directory to mesh file.
    """

    mesh, mf_boundaries, association_table = import_mesh(prefix=mesh_name, subdomains=False, directory=mesh_directory)

    X = mesh.coordinates()
    X = [X[:, i] for i in range(mesh.geometry().dim())]

    plt.figure(figsize=(20, 10))
    plt.triplot(mesh2string(mesh))
    plt.title('mesh at {}'.format(mesh_name))
    plt.axis('off')
    plt.show()


if __name__ == "__main__":
    simulation_type = 'ellipse'
    simulation_res = 'has'
    index = 447
    time_stamp = 719
    mesh_file = '{}_{}'.format(simulation_type, index)
    mesh_directory = "mesh/{}/{}/".format(simulation_type, simulation_res)
    solution_directory = 'solution/{}_{}_{}/data/1/Timeseries'.format(simulation_type, index, simulation_res)
    t = 800

    # plot_solution(mesh_file, mesh_directory, solution_directory, t)
    plot_prediction(mesh_file, mesh_directory, 'prediction/prediction_result_{}_{}_{}_model_675.npz'.format(simulation_type, index, time_stamp))
    #                                                                                                                                                                                                plot_prediction(mesh_file, mesh_directory, 'prediction/prediction_result_{}_{}_{}_model_175.npz'.format(simulation_type, index, time_stamp))
    # plot_mesh(mesh_file, mesh_directory)