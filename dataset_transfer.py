from dolfin import *
from fenics import *
from utils import *

if __name__ == "__main__":
    for i in range(16, 17):
        mesh_name = 'circle_{}'.format(i)
        type = 'circle'
        # check if calculation is stable
        # if ensure_stable_calculation(mesh_name, type):
        #     print('Calculation {} is stable.'.format(mesh_name))
            # read and extract data
        read_timeseries_to_npy(mesh_name, type)
        # else:
        #     print('Calculation {} is unstable.'.format(mesh_name))

    # for i in range(1000):
    #     mesh_name = 'ellipse_{}'.format(i)
    #     type = 'ellipse'
    #     # check if calculation is stable
    #     if ensure_stable_calculation(mesh_name, type):
    #         print('Calculation {} is stable.'.format(mesh_name))
    #         # read and extract data
    #         read_timeseries_to_npy(mesh_name, type)
    #     else:
    #         print('Calculation {} is unstable.'.format(mesh_name))