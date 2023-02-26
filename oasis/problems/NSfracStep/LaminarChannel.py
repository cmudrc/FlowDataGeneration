from __future__ import print_function
__author__ = "Mikael Mortensen <mikaem@math.uio.no>"
__date__ = "2013-06-25"
__copyright__ = "Copyright (C) 2013 " + __author__
__license__ = "GNU Lesser GPL version 3 or any later version"

from ..NSfracStep import *
import numpy as np
from numpy import pi, arctan, array
from utils import import_mesh
import global_variables
set_log_active(False)

prefix = global_variables.get_value('meshname')
directory = global_variables.get_value('meshdir')

# Override some problem specific parameters
def problem_parameters(NS_parameters, NS_expressions, **NS_namespace):
    nu = 0.001
    Re = 1. / nu
    Um = 4
    geo_spec = np.load(directory + '/' + prefix + '.npz')
    L = geo_spec['domain_size'][0]
    H = geo_spec['domain_size'][1]
    NS_parameters.update(dict(
        nu=nu,
        L=L,
        H=H,
        Re=Re,
        Um=Um,
        # Nx=40,
        # Ny=40,
        max_iter=1000,
        velocity_degree=2,
        use_krylov_solvers=True,
        krylov_solvers=dict(monitor_convergence=True)))

    # NS_expressions.update(dict(constrained_domain=PeriodicDomain(L)))


# Load mesh
def mesh(**params):
    mesh, mf_boundaries, association_table = import_mesh(prefix=prefix, subdomains=False, directory=directory)
    return mesh


def create_bcs(V, Q, L, H, Um, sys_comp, **NS_namespace):
    # def walls(x, on_boundary):
    #     return (on_boundary and (x[0] > 1e-8 and x[0] < L - 1e-8))
    inlet = Expression(
        "4.*{0}*x[1]*({1}-x[1])/pow({1}, 2)".format(Um, H), degree=2)
    Walls = AutoSubDomain(lambda x, on_bnd: on_bnd and x[0] > 1e-8 and x[0] < L - 1e-8)
    Inlet = AutoSubDomain(lambda x, on_bnd: on_bnd and x[0] < 1e-8)
    Outlet = AutoSubDomain(lambda x, on_bnd: on_bnd and x[0] > L - 1e-8)

    bcs = dict((ui, []) for ui in sys_comp)
    bc00 = DirichletBC(V, inlet, Inlet)
    bc01 = DirichletBC(V, 0, Inlet)
    bc1 = DirichletBC(V, 0, Walls)
    bcp = DirichletBC(Q, 0, Outlet)
    return dict(u0=[bc00, bc1],
                u1=[bc01, bc1],
                p=[bcp])


def body_force(Re, **NS_namespace):
    return Constant((2. / Re, 0.))


# def reference(Re, t, num_terms=100):
#     u = 1.0
#     c = 1.0
#     for n in range(1, 2 * num_terms, 2):
#         a = 32. / (pi**3 * n**3)
#         b = (0.25 / Re) * pi**2 * n**2
#         c = -c
#         u += a * exp(-b * t) * c
#     return u


def temporal_hook(tstep, q_, t, Re, L, **NS_namespace):
    if tstep % 20 == 0:
        plot(q_['u0'])
    try:
        # point is found on one processor, the others pass
        u_computed = q_['u0'](array([L, 0.]))
    except:
        pass
