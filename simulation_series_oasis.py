import os
from fenics import *
from dolfin import *

# from mshr import *


if __name__ == "__main__":
    T = 10            # final time
    num_steps = 1000   # number of time steps
    dt = T / num_steps # time step size
    mu = 0.001         # dynamic viscosity
    rho = 1            # density
    
    directory = 'mesh/ellipse/'
    for i in range(1000):
        mesh_file_name = 'ellipse_{}'.format(i)
        mesh_directory = os.path.join(directory, 'las')
        save_file_dir = 'solution/ellipse_{}_las'.format(i)
        if not os.path.exists(save_file_dir):
            os.makedirs(save_file_dir)

        os.system('python NSfracStep.py problem=Cylinder meshname={} meshdir={} T={} dt={} output_timeseries_as_vector=False folder={}'.format(mesh_file_name, mesh_directory, T, dt, save_file_dir))
        # os.system('python NSfracStep.py problem=LaminarChannel meshname={} meshdir={} T={} dt={} output_timeseries_as_vector=False folder={}'.format(mesh_file_name, mesh_directory, T, dt, save_file_dir))