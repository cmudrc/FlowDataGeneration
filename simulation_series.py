import os
from ast import Interactive
from fenics import *
from dolfin import *
# from mshr import *
import matplotlib.pyplot as plt
from utils import *


# Define symmetric gradient
def epsilon(u):
    return sym(nabla_grad(u))

# Define stress tensor
def sigma(u, p):
    return 2*mu*epsilon(u) - p*Identity(len(u))

def solve_geometry_domain(mesh_file_name, directory, save_file_dir, inflow_profile, num_steps, dt, mu, rho):
    '''
    Solves target geometry domain with FEniCS.
    mesh_file_name: the prefix of mesh file, for 'cylinder.xdmf' simply input 'cylinder'. DO NOT INCLUDE directory here.
    directory: dir to mesh file
    save_file_dir: dir to output file. Output would be stored by nodal values in 'velocity_series' and 'pressure_series' separately.
    inflow_profile: inflow profile defined in FEniCS format, or a string defining the inflow function, e.g. '4.0*1.5*x[1]*(10 - x[1]) / 100'  
    '''

    mesh, mf_boundaries, association_table = import_mesh(prefix=mesh_file_name, subdomains=False, directory=directory)

    inflow_id = association_table['inflow']
    outflow_id = association_table['outflow']
    walls_id = association_table['walls']
    cylinder_id = association_table['circle']

    # Define function spaces
    V = VectorFunctionSpace(mesh, 'P', 2)
    Q = FunctionSpace(mesh, 'P', 1)

    # Define inflow profile
    inflow_profile = inflow_profile

    # Define boundary conditions
    bcu_inflow = DirichletBC(V, Expression(inflow_profile, degree=2), mf_boundaries, inflow_id)
    # bcu_inflow = DirichletBC(V, Expression(inflow_profile, degree=2), df_inflow)
    bcu_walls = DirichletBC(V, Constant((0, 0)), mf_boundaries, walls_id)
    bcu_cylinder = DirichletBC(V, Constant((0, 0)), mf_boundaries, cylinder_id)
    bcp_outflow = DirichletBC(Q, Constant(0), mf_boundaries, outflow_id)
    bcu = [bcu_inflow, bcu_walls, bcu_cylinder]
    bcp = [bcp_outflow]

    # Define trial and test functions
    u = TrialFunction(V)
    v = TestFunction(V)
    p = TrialFunction(Q)
    q = TestFunction(Q)

    # Define functions for solutions at previous and current time steps
    u_n = Function(V)
    u_  = Function(V)
    p_n = Function(Q)
    p_  = Function(Q)

    # Define expressions used in variational forms
    U = 0.5*(u_n + u)
    n = FacetNormal(mesh)
    f = Constant((0, 0))
    k = Constant(dt)
    mu = Constant(mu)
    rho = Constant(rho)

    # Define variational problem for step 1
    F1 = rho*dot((u - u_n) / k, v)*dx \
    + rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx \
    + inner(sigma(U, p_n), epsilon(v))*dx \
    + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds \
    - dot(f, v)*dx
    a1 = lhs(F1)
    L1 = rhs(F1)

    # Define variational problem for step 2
    a2 = dot(nabla_grad(p), nabla_grad(q))*dx
    L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

    # Define variational problem for step 3
    a3 = dot(u, v)*dx
    L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx

    # Assemble matrices
    A1 = assemble(a1)
    A2 = assemble(a2)
    A3 = assemble(a3)

    # Apply boundary conditions to matrices
    [bc.apply(A1) for bc in bcu]
    [bc.apply(A2) for bc in bcp]

    # Create time series (for use in reaction_system.py)
    timeseries_u = TimeSeries(os.path.join(save_file_dir, 'velocity_series'))
    timeseries_p = TimeSeries(os.path.join(save_file_dir, 'pressure_series'))

    # Create progress bar
    progress = Progress('Time-stepping')
    # set_log_level(PROGRESS)

    # Time-stepping
    t = 0
    for n in range(num_steps):

        # Update current time
        t += dt

        # Step 1: Tentative velocity step
        b1 = assemble(L1)
        [bc.apply(b1) for bc in bcu]
        solve(A1, u_.vector(), b1, 'bicgstab', 'hypre_amg')

        # Step 2: Pressure correction step
        b2 = assemble(L2)
        [bc.apply(b2) for bc in bcp]
        solve(A2, p_.vector(), b2, 'bicgstab', 'hypre_amg')

        # Step 3: Velocity correction step
        b3 = assemble(L3)
        solve(A3, u_.vector(), b3, 'cg', 'sor')

        # Save nodal values to file
        timeseries_u.store(u_.vector(), t)
        timeseries_p.store(p_.vector(), t)

        # Update previous solution
        u_n.assign(u_)
        p_n.assign(p_)

        # Update progress bar
        # progress.update(t / T)
        print('u max:', u_.vector().get_local().max())

if __name__ == "__main__":
    T = 10            # final time
    num_steps = 8000   # number of time steps
    dt = T / num_steps # time step size
    mu = 0.001         # dynamic viscosity
    rho = 1            # density
    
    directory = 'io_operations/data221120/'
    inflow_profile = ('4.0*1.5*x[1]*(10 - x[1]) / 100', '0')
    for i in range(2, 3):
        mesh_file_name = 'circle_{}'.format(i)
        mesh_directory = os.path.join(directory, 'has')
        save_file_dir = 'navier_stokes_cylinder/sol221111_{}_has'.format(i)
        if not os.path.exists(save_file_dir):
            os.makedirs(save_file_dir)

        solve_geometry_domain(mesh_file_name, mesh_directory, save_file_dir, inflow_profile, num_steps, dt, mu, rho)