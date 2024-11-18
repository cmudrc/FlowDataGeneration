__author__ = "Mikael Mortensen <mikaem@math.uio.no>"
__date__ = "2014-04-10"
__copyright__ = "Copyright (C) 2014 " + __author__
__license__ = "GNU Lesser GPL version 3 or any later version"

from dolfin import Mesh, AutoSubDomain, near
from utils import import_mesh
from oasis.problems.NSfracStep import *
import global_variables
import os
import platform
import numpy as np

# if not os.path.isfile("cylinder.xml"):
#     if platform.system() == "Linux":
#         os.system("wget -O cylinder.xml https://www.dropbox.com/s/d78g4cyjxl3ylay/cylinder.xml?dl=0")
#     elif platform.system() == "Darwin":
#         os.system("curl -L https://www.dropbox.com/s/d78g4cyjxl3ylay/cylinder.xml?dl=0 -o cylinder.xml")
#     else:
#         raise ImportError("Could not determine platform")

    # try:
        #os.system("gmsh mesh/cylinder.geo -2 -o mesh/cylinder.msh")
        #os.system("dolfin-convert mesh/cylinder.msh mesh/cylinder.xml")
        #os.system("rm mesh/cylinder.msh")
    # except RuntimeError:
        #os.system("wget -O cylinder.xml https://www.dropbox.com/s/d78g4cyjxl3ylay/cylinder.xml?dl=0")
        ##raise "Gmsh is required to run this demo"

# mesh = Mesh("cylinder.xml")

prefix = global_variables.get_value('meshname')
directory = global_variables.get_value('meshdir')
type = prefix.split('_')[0]

mesh, mf_boundaries, association_table = import_mesh(prefix=prefix, subdomains=False, directory=directory)
geo_spec = np.load(directory + '/' + prefix + '.npz')

H = geo_spec['domain_size'][1]
L = geo_spec['domain_size'][0]

D = 2 * np.sum(geo_spec['size'])

center_x = geo_spec['position'][0]
center_y = geo_spec['position'][1]
cases = {
    1: {'Um': 500,
        'Re': 500 * D / 0.001},

    2: {'Um': 10,
        'Re': 1000.0}
}

# Specify boundary conditions

inflow_id = association_table['inflow']
outflow_id = association_table['outflow']
# walls_id = association_table['wall']
cylinder_id = association_table[type]

# Inlet = inflow_id
# Wall = walls_id
# Cyl = cylinder_id
# Outlet = outflow_id

Inlet = AutoSubDomain(lambda x, on_bnd: on_bnd and x[0] < 1e-8)
Wall = AutoSubDomain(lambda x, on_bnd: on_bnd and near(x[1] * (H - x[1]), 0))
Cyl = AutoSubDomain(lambda x, on_bnd: (on_bnd and x[0] > 1e-4 and x[0] < H-1e-4
                                      and x[1] < H-1e-4 and x[1] > 1e-4))
Outlet = AutoSubDomain(lambda x, on_bnd: on_bnd and x[0] > L - 1e-8)


# Overload post_import_problem to choose between the two cases
def post_import_problem(NS_parameters, commandline_kwargs, **NS_namespace):
    """ Choose case - case could be defined through command line."""
    NS_parameters.update(commandline_kwargs)
    case = NS_parameters['case'] if 'case' in NS_parameters else 1
    Um = cases[case]["Um"]
    Re = cases[case]["Re"]
    Umean = 2. / 3. * Um
    nu = Umean * D / Re
    NS_parameters.update(nu=nu, Re=Re, Um=Um, Umean=Umean)

    return NS_parameters
