import os
from fenics import *
from dolfin import *

# from mshr import *


if __name__ == "__main__":
    T = 10            # final time
    num_steps = 10000   # number of time steps
    dt = T / num_steps # time step size
    mu = 0.001         # dynamic viscosity
    rho = 1            # density
    
    directory = 'io_operations/data221206/'
    for i in range(2, 3):
        mesh_file_name = 'channel_{}'.format(i)
        mesh_directory = os.path.join(directory, 'has')
        save_file_dir = 'navier_stokes_cylinder/sol221129_{}_has'.format(i)
        if not os.path.exists(save_file_dir):
            os.makedirs(save_file_dir)

        # os.system('python NSfracStep.py problem=Cylinder meshname={} meshdir={} T={} dt={} output_timeseries_as_vector=False folder={}'.format(mesh_file_name, mesh_directory, T, dt, save_file_dir))
        os.system('python NSfracStep.py problem=LaminarChannel meshname={} meshdir={} T={} dt={} output_timeseries_as_vector=False folder={}'.format(mesh_file_name, mesh_directory, T, dt, save_file_dir))