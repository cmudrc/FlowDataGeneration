import meshio
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import matplotlib.ticker as ticker

from dolfin import *
from fenics import *
from fenicstools.Interpolation import interpolate_nonmatching_mesh
from configparser import ConfigParser
try:
    from dolfin import XDMFFile, Mesh, MeshValueCollection
    from dolfin.cpp.mesh import MeshFunctionSizet
except ImportError:
    print("Could not import dolfin. Continuing without Dolfin support.")


def msh2xdmf(mesh_name, dim=2, directory="."):
    """
    Function converting a MSH mesh into XDMF files.
    The XDMF files are:
        - "domain.xdmf": the domain;
        - "boundaries.xdmf": the boundaries physical groups from GMSH;
    """

    # Get the mesh name has prefix
    prefix = mesh_name.split('.')[0]
    # Read the input mesh
    msh = meshio.read("{}/{}".format(directory, mesh_name))
    # Generate the domain XDMF file
    export_domain(msh, dim, directory, prefix)
    # Generate the boundaries XDMF file
    export_boundaries(msh, dim, directory, prefix)
    # Export association table
    export_association_table(msh, prefix, directory)


def export_domain(msh, dim, directory, prefix):
    """
    Export the domain XDMF file as well as the subdomains values.
    """
    # Set cell type
    if dim == 2:
        cell_type = "triangle"
    elif dim == 3:
        cell_type = "tetra"
    # Generate the cell block for the domain cells
    for i in msh.cells:
        if i.type == cell_type:
            data_array = i.data 
    # data_array = [arr for (t, arr) in msh.cells if t == cell_type]
    
    if len(data_array) == 0:
        print("WARNING: No domain physical group found.")
        return
    else:
        # data = np.concatenate(data_array) # Use this expression if more than 1 domain
        data = data_array
    cells = [
        meshio.CellBlock(
            cell_type=cell_type,
            data=data,
        )
    ]
    # Generate the domain cells data (for the subdomains)
    try:
        cell_data = {
            "subdomains": [
                np.concatenate(
                    [
                        msh.cell_data["gmsh:physical"][i]
                        for i, cellBlock in enumerate(msh.cells)
                        if cellBlock.type == cell_type
                    ]
                )
            ]
        }
    except KeyError:
        raise ValueError(
            """
            No physical group found for the domain.
            Define the domain physical group.
                - if dim=2, the domain is a surface
                - if dim=3, the domain is a volume
            """
        )

    # Generate a meshio Mesh for the domain
    domain = meshio.Mesh(
        points=msh.points[:, :dim],
        cells=cells,
        cell_data=cell_data,
    )
    # Export the XDMF mesh of the domain
    meshio.write(
        "{}/{}_{}".format(directory, prefix, "domain.xdmf"),
        domain,
        file_format="xdmf"
    )


def export_boundaries(msh, dim, directory, prefix):
    """
    Export the boundaries XDMF file.
    """
    # Set the cell type
    if dim == 2:
        cell_type = "line"
    elif dim == 3:
        cell_type = "triangle"
    # Generate the cell block for the boundaries cells
    # data_array = [arr for (t, arr) in msh.cells if t == cell_type]
    data_array = []
    for i in msh.cells:
        if i.type == cell_type:
            data_array.append(i.data) 
    if len(data_array) == 0:
        print("WARNING: No boundary physical group found.")
        return
    else:
        data = np.concatenate(data_array)
        # data = data_array
    boundaries_cells = [
        meshio.CellBlock(
            cell_type=cell_type,
            data=data,
        )
    ]
    # Generate the boundaries cells data
    cell_data = {
        "boundaries": [
            np.concatenate(
                [
                    msh.cell_data["gmsh:physical"][i]
                    for i, cellBlock in enumerate(msh.cells)
                    if cellBlock.type == cell_type
                ]
            )
        ]
    }
    # Generate the meshio Mesh for the boundaries physical groups
    boundaries = meshio.Mesh(
        points=msh.points[:, :dim],
        cells=boundaries_cells,
        cell_data=cell_data,
    )
    # Export the XDMF mesh of the lines boundaries
    meshio.write(
        "{}/{}_{}".format(directory, prefix, "boundaries.xdmf"),
        boundaries,
        file_format="xdmf"
    )


def export_association_table(msh, prefix='mesh', directory='.', verbose=True):
    """
    Display the association between the physical group label and the mesh
    value.
    """
    # Create association table
    association_table = {}

    # Display the correspondance
    formatter = "|{:^20}|{:^20}|"
    topbot = "+{:-^41}+".format("")
    separator = "+{:-^20}+{:-^20}+".format("", "")

    # Display
    if verbose:
        print('\n' + topbot)
        print(formatter.format("GMSH label", "MeshFunction value"))
        print(separator)

    for label, arrays in msh.cell_sets.items():
        # Get the index of the array in arrays
        for i, array in enumerate(arrays):
            if array.size != 0:
                index = i
        # Added check to make sure that the association table
        # doesn't try to import irrelevant information.
        if label != "gmsh:bounding_entities":
            value = msh.cell_data["gmsh:physical"][index][0]
            # Store the association table in a dictionnary
            association_table[label] = value
            # Display the association
            if verbose:
                print(formatter.format(label, value))
    if verbose:
        print(topbot)
    # Export the association table
    file_content = ConfigParser()
    file_content["ASSOCIATION TABLE"] = association_table
    file_name = "{}/{}_{}".format(directory, prefix, "association_table.ini")
    with open(file_name, 'w') as f:
        file_content.write(f)


def import_mesh(
        prefix="mesh",
        subdomains=False,
        dim=2,
        directory=".",
):
    """Function importing a dolfin mesh.

    Arguments:
        prefix (str, optional): mesh files prefix (eg. my_mesh.msh,
            my_mesh_domain.xdmf, my_mesh_bondaries.xdmf). Defaults to "mesh".
        subdomains (bool, optional): True if there are subdomains. Defaults to
            False.
        dim (int, optional): dimension of the domain. Defaults to 2.
        directory (str, optional): directory of the mesh files. Defaults to ".".

    Output:
        - dolfin Mesh object containing the domain;
        - dolfin MeshFunction object containing the physical lines (dim=2) or
            surfaces (dim=3) defined in the msh file and the sub-domains;
        - association table
    """
    # Set the file name
    domain = "{}_domain.xdmf".format(prefix)
    boundaries = "{}_boundaries.xdmf".format(prefix)

    # create 2 xdmf files if not converted before
    if not os.path.exists("{}/{}".format(directory, domain)) or \
       not os.path.exists("{}/{}".format(directory, boundaries)):
        msh2xdmf("{}.msh".format(prefix), dim=dim, directory=directory)

    # Import the converted domain
    mesh = Mesh()
    with XDMFFile("{}/{}".format(directory, domain)) as infile:
        infile.read(mesh)
    # Import the boundaries
    boundaries_mvc = MeshValueCollection("size_t", mesh, dim=dim)
    with XDMFFile("{}/{}".format(directory, boundaries)) as infile:
        infile.read(boundaries_mvc, 'boundaries')
    boundaries_mf = MeshFunctionSizet(mesh, boundaries_mvc)
    # Import the subdomains
    if subdomains:
        subdomains_mvc = MeshValueCollection("size_t", mesh, dim=dim)
        with XDMFFile("{}/{}".format(directory, domain)) as infile:
            infile.read(subdomains_mvc, 'subdomains')
        subdomains_mf = MeshFunctionSizet(mesh, subdomains_mvc)
    # Import the association table
    association_table_name = "{}/{}_{}".format(
        directory, prefix, "association_table.ini")
    file_content = ConfigParser()
    file_content.read(association_table_name)
    association_table = dict(file_content["ASSOCIATION TABLE"])
    # Convert the value from string to int
    for key, value in association_table.items():
        association_table[key] = int(value)
    # Return the Mesh and the MeshFunction objects
    if not subdomains:
        return mesh, boundaries_mf, association_table
    else:
        return mesh, boundaries_mf, subdomains_mf, association_table


def read_timeseries_to_npy(mesh_name, type):
    """
    Read Oasis TimeSeries and export to numpy arrays.
    Input:
        mesh_name: mesh prefix to read
        type: 'circle' or 'ellipses' or 'channel'
    Output:
        None, value of all timesteps stored to numpy files
    """
    mesh_l, mf_boundaries_l, association_table_l = import_mesh(prefix=mesh_name, directory='mesh/{}/las'.format(type))
    mesh_h, mf_boundaries_h, association_table_h = import_mesh(prefix=mesh_name, directory='mesh/{}/has'.format(type))

    geometry_index = mesh_name

    gdim = mesh_l.geometry().dim()
    tdim = mesh_l.topology().dim()

    V = FunctionSpace(mesh_h, 'CG', 2)
    Q = FunctionSpace(mesh_h, 'CG', 1)
    u0_ = Function(V)
    u1_ = Function(V)
    p_ = Function(Q)

    Vl = FunctionSpace(mesh_l, 'CG', 2)
    Ql = FunctionSpace(mesh_l, 'CG', 1)
    u0_l = Function(Vl)
    u1_l = Function(Vl)
    p_l = Function(Ql)

    # nv = mesh_l.num_vertices()
    X = mesh_h.coordinates()
    X = [X[:, i] for i in range(gdim)]

    # Store mesh edges
    lines = np.zeros((2*mesh_h.num_edges(), 2))
    line_length = np.zeros(2*mesh_h.num_edges())
    cells = np.array(mesh_h.cells())
    for i, edge in enumerate(edges(mesh_h)):
        lines[2*i, :] = edge.entities(0)
        lines[2*i+1, :] = np.flipud(edge.entities(0))
        line_length[2*i] = edge.length()
        line_length[2*i+1] = edge.length()

    # Read solution
    x = X[0]
    y = X[1]

    velocity_x = TimeSeries('solution/{}_has/data/1/Timeseries/u0_from_tstep_0'.format(mesh_name))
    velocity_y = TimeSeries('solution/{}_has/data/1/Timeseries/u1_from_tstep_0'.format(mesh_name))
    pressure = TimeSeries('solution/{}_has/data/1/Timeseries/p_from_tstep_0'.format(mesh_name))

    velocity_xl = TimeSeries('solution/{}_las/data/1/Timeseries/u0_from_tstep_0'.format(mesh_name))
    velocity_yl = TimeSeries('solution/{}_las/data/1/Timeseries/u1_from_tstep_0'.format(mesh_name))
    pressure_l = TimeSeries('solution/{}_las/data/1/Timeseries/p_from_tstep_0'.format(mesh_name))
    for t in range(1000, 1001, 1):
        time_index = t
        velocity_x.retrieve(u0_.vector(), t)
        velocity_y.retrieve(u1_.vector(), t)
        pressure.retrieve(p_.vector(), t)
        # # plot velocity contour  
        # plot(u0_, title='u0')
        # plt.savefig('data/has/{}_has_{}_u0.png'.format(mesh_name, t))
        # plt.figure()
        # plot(u1_, title='u1')
        # plt.savefig('data/has/{}_has_{}_u1.png'.format(mesh_name, t))
        # plt.show()
        w0 = u0_.compute_vertex_values(mesh_h)
        w1 = u1_.compute_vertex_values(mesh_h)
        C = p_.compute_vertex_values(mesh_h)

        result_hx = w0
        result_hy = w1
        result_hp = C

        if os.path.exists('data/has'):
            np.savez('data/has/{}_has_{}'.format(mesh_name, t), ux=result_hx, uy=result_hy, p=result_hp)
        else:
            os.makedirs('data/has')
            np.savez('data/has/{}_has_{}'.format(mesh_name, t), ux=result_hx, uy=result_hy, p=result_hp)

        velocity_xl.retrieve(u0_l.vector(), t)
        velocity_yl.retrieve(u1_l.vector(), t)
        pressure_l.retrieve(p_l.vector(), t)
        
        ###############################################
        # plot velocity contour of las and has and overlay with mesh
        fig, ax = plt.subplots(2, 1, figsize=(18, 12))
        plt.rcParams.update({'font.size': 20})
        w = np.sqrt(w0**2 + w1**2)
        ax[1].tricontourf(x, y, cells, w, cmap='jet', levels=100)
        tri_h = Triangulation(x, y, cells)
        # ax[1].triplot(tri_h, 'k-', linewidth=0.5)
        ax[1].axis('off')
        ax[1].set_title('(b) High resolution simulation', y=-0.1)

        w0_l = u0_l.compute_vertex_values(mesh_l)
        w1_l = u1_l.compute_vertex_values(mesh_l)
        w_l = np.sqrt(w0_l**2 + w1_l**2)
        x_l = mesh_l.coordinates()[:, 0]
        y_l = mesh_l.coordinates()[:, 1]
        cells_l = np.array(mesh_l.cells())
        ax[0].tricontourf(x_l, y_l, cells_l, w_l, cmap='jet', levels=100)
        tri_l = Triangulation(x_l, y_l, cells_l)
        ax[0].triplot(tri_l, 'k-', linewidth=0.5)
        ax[0].axis('off')
        ax[0].set_title('(a) Low resolution simulation', y=-0.1)
        fig.tight_layout()
        cb1 = fig.colorbar(ax[1].tricontourf(x, y, cells, w, cmap='jet', levels=100), ax=ax.ravel().tolist(), shrink=0.95, location='right', pad=0.05, label='Velocity magnitude (m/s)')
        cb1.locator = ticker.MaxNLocator(nbins=6)
        # cb1.update_ticks()
        cb1.set_ticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
        ax[1].triplot(tri_h, 'k-', linewidth=0.5)
        plt.savefig('data/{}_compare_{}_velocity.png'.format(mesh_name, t))


        # plot pressure contour of las and has and overlay with mesh
        fig, ax = plt.subplots(2, 1, figsize=(18, 12))
        plt.rcParams.update({'font.size': 20})
        ax[1].tricontourf(x, y, cells, C, cmap='jet', levels=100)
        tri_h = Triangulation(x, y, cells)
        
        ax[1].axis('off')
        ax[1].set_title('(b) High resolution simulation', y=-0.1)

        C_l = p_l.compute_vertex_values(mesh_l)
        ax[0].tricontourf(x_l, y_l, cells_l, C_l, cmap='jet', levels=100)
        tri_l = Triangulation(x_l, y_l, cells_l)
        ax[0].triplot(tri_l, 'k-', linewidth=0.5)
        ax[0].axis('off')
        ax[0].set_title('(a) Low resolution simulation', y=-0.1)
        fig.tight_layout()
        cb2 = fig.colorbar(ax[1].tricontourf(x, y, cells, C, cmap='jet', levels=100), ax=ax.ravel().tolist(), shrink=0.95, location='right', pad=0.05, label='Pressure (Pa)')
        cb2.locator = ticker.MaxNLocator(nbins=8)
        # cb2.update_ticks()
        cb2.set_ticks([-4.5, -3.0, -1.0, 1.0, 3.0, 5.0, 7.0, 9.0])
        ax[1].triplot(tri_h, 'k-', linewidth=0.5)
        plt.savefig('data/{}_compare_{}_pressure.png'.format(mesh_name, t))

        ###############################################

        u0_l_ = interpolate_nonmatching_mesh(u0_l, V)
        u1_l_ = interpolate_nonmatching_mesh(u1_l, V)
        p_l_ = interpolate_nonmatching_mesh(p_l, Q)
                
        # plt.figure()
        # plot(u0_l_, title='u0')
        # plt.savefig('data/las/{}_las_{}_u0_interpolate.png'.format(mesh_name, t))
        # plt.figure()
        # plot(u1_l_, title='u1')
        # plt.savefig('data/las/{}_las_{}_u1_interpolate.png'.format(mesh_name, t))

        u0_l_ = u0_l_.compute_vertex_values(mesh_h)
        u1_l_ = u1_l_.compute_vertex_values(mesh_h)
        p_l_ = p_l_.compute_vertex_values(mesh_h)

        result_lx = u0_l_
        result_ly = u1_l_
        result_lp = p_l_

        if os.path.exists('data/las'):
            np.savez('data/las/{}_las_{}'.format(mesh_name, t), ux=result_lx, uy=result_ly, p=result_lp, geometry_index=geometry_index, time_index=time_index)
        else:
            os.makedirs('data/las')
            np.savez('data/las/{}_las_{}'.format(mesh_name, t), ux=result_lx, uy=result_ly, p=result_lp, geometry_index=geometry_index, time_index=time_index)

    # save mesh
    if os.path.exists('data/mesh'):
        np.savez('data/mesh/{}'.format(mesh_name), x=x, y=y, edges=lines, edge_properties=line_length, cells=cells, geometry_index=geometry_index)
    else:
        os.makedirs('data/mesh')
        np.savez('data/mesh/{}'.format(mesh_name), x=x, y=y, edges=lines, edge_properties=line_length, cells=cells, geometry_index=geometry_index)
    

def ensure_stable_calculation(mesh_name, type):
    """
    A small tool to check if the calculation is stable by examining the value at 1000th timestep.
    Input:
        mesh_name: mesh prefix to read
        type: 'circle' or 'ellipses' or 'channel'
    Output:
        flag: True if stable, False if not
    """
    mesh_l, mf_boundaries_l, association_table_l = import_mesh(prefix=mesh_name, directory='mesh/{}/las'.format(type))
    mesh_h, mf_boundaries_h, association_table_h = import_mesh(prefix=mesh_name, directory='mesh/{}/has'.format(type))

    gdim = mesh_l.geometry().dim()
    tdim = mesh_l.topology().dim()

    Q = FunctionSpace(mesh_h, 'CG', 1)
    p_ = Function(Q)
    p1_ = Function(Q)

    Ql = FunctionSpace(mesh_l, 'CG', 1)
    pl_ = Function(Ql)
    p1l_ = Function(Ql)

    X = mesh_l.coordinates()
    X = [X[:, i] for i in range(gdim)]

    pressure = TimeSeries('solution/{}_has/data/1/Timeseries/p_from_tstep_0'.format(mesh_name))
    pressure_l = TimeSeries('solution/{}_las/data/1/Timeseries/p_from_tstep_0'.format(mesh_name))

    # check if 1000 timestep can be retrieved
    
    pressure.retrieve(p_.vector(), 998)
    pressure.retrieve(p1_.vector(), 999)
    pressure_l.retrieve(pl_.vector(), 998)
    pressure_l.retrieve(p1l_.vector(), 999)
    if np.allclose(p_.compute_vertex_values(mesh_h), p1_.compute_vertex_values(mesh_h)) or np.allclose(pl_.compute_vertex_values(mesh_l), p1l_.compute_vertex_values(mesh_l)):
        flag = False
    else:
        flag = True
    
    return flag
    

def read_mesh_to_npy(mesh_name, type):
    """
    Read mesh and solution from .xdmf file and save them as .npy file
    Input:
        mesh_name: mesh prefix to read
        type: 'circle' or 'ellipses' or 'channel'
    Output:
        None
    """
    mesh_l, mf_boundaries_l, association_table_l = import_mesh(prefix=mesh_name, directory='mesh/{}/las'.format(type))
    mesh_h, mf_boundaries_h, association_table_h = import_mesh(prefix=mesh_name, directory='mesh/{}/has'.format(type))

    gdim = mesh_l.geometry().dim()

    X = mesh_l.coordinates()
    X = [X[:, i] for i in range(gdim)]

    # Store mesh edges
    lines = np.zeros((2*mesh_l.num_edges(), 2))
    line_length = np.zeros(2*mesh_l.num_edges())
    for i, edge in enumerate(edges(mesh_l)):
        lines[2*i, :] = edge.entities(0)
        lines[2*i+1, :] = np.flipud(edge.entities(0))
        line_length[2*i] = edge.length()
        line_length[2*i+1] = edge.length()

    # Read solution
    x = X[0]
    y = X[1]

    # save mesh
    if os.path.exists('data/mesh/las'):
        np.savez('data/mesh/las/{}'.format(mesh_name), x=x, y=y, edges=lines, edge_properties=line_length)
    else:
        os.makedirs('data/mesh/las', exist_ok=True)
        np.savez('data/mesh/las/{}'.format(mesh_name), x=x, y=y, edges=lines, edge_properties=line_length)
   
    X = mesh_h.coordinates()
    X = [X[:, i] for i in range(gdim)]

    # Store mesh edges
    lines = np.zeros((2*mesh_h.num_edges(), 2))
    line_length = np.zeros(2*mesh_h.num_edges())
    for i, edge in enumerate(edges(mesh_h)):
        lines[2*i, :] = edge.entities(0)
        lines[2*i+1, :] = np.flipud(edge.entities(0))
        line_length[2*i] = edge.length()
        line_length[2*i+1] = edge.length()

    # Read solution
    x = X[0]
    y = X[1]

    # save mesh
    if os.path.exists('data/mesh/has'):
        np.savez('data/mesh/has/{}'.format(mesh_name), x=x, y=y, edges=lines, edge_properties=line_length)
    else:
        os.makedirs('data/mesh/has', exist_ok=True)
        np.savez('data/mesh/has/{}'.format(mesh_name), x=x, y=y, edges=lines, edge_properties=line_length)
   