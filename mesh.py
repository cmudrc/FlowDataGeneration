
import os
from cmath import pi
import numpy as np
import gmsh
from utils import *
import random

def createFreeFlowMesh(domain_size, type, position, size, mesh_size = 0.2, filename = '', directory = ''):
    '''
    Generate mesh
    Input:
        domain_size: Absolute size in scale (relative to mesh size)
        type: cylinder, ellipse, airfoil
        position: x and y coordinates of shape center. For type = ellipse add the angle of attack as the third element in this array. 
        size: Passed in accordance with type. For circle objects it should indicate the radius of the circle. For ellipses pass the long and short axis 
            in form of numpy array (e.g. np.array([5, 3])), note that the first element (long axis) must be larger than the second element. 
        mesh_size: Mesh size
        filename: Name of the mesh file
        directory: Directory where the mesh file is saved
    Output:
        None
    '''
    gmsh.initialize()
    gmsh.clear()
    gmsh.model.add("t3")

    # Initialize simulation domain
    gmsh.model.occ.addPoint(domain_size[0], 0, 0, tag=1)
    gmsh.model.occ.addPoint(domain_size[0], domain_size[1], 0, tag=2)
    gmsh.model.occ.addPoint(0, 0, 0, tag=3)
    gmsh.model.occ.addPoint(0, domain_size[1], 0, tag=4)
 
    gmsh.model.occ.addLine(3, 1, 1)
    gmsh.model.occ.addLine(3, 4, 2)
    gmsh.model.occ.addLine(1, 2, 3)
    gmsh.model.occ.addLine(2, 4, 4)
    gmsh.model.occ.addCurveLoop([1, 2, 3, 4], tag=1)

    # gmsh.model.occ.addCircle(5, 6, 7, tag=5)
    if(type == 'circle'):
        gmsh.model.occ.addPoint(position[0]+size, position[1], 0, tag=5) # Control mesh size on boundaries
        gmsh.model.occ.addCircle(position[0], position[1], 0, size, tag=5, angle1=0, angle2=2*pi)
        gmsh.model.occ.addCurveLoop([5], tag=2)
        gmsh.model.occ.addPlaneSurface([1, 2], tag=1)
        gmsh.model.occ.mesh.setSize(gmsh.model.occ.getEntities(0), mesh_size)

    elif(type == 'ellipse'):
        xAxis = np.array([1, np.arctan(position[2] * np.pi / 180), 0])
        zAxis = np.array([0, 0, 1])
        gmsh.model.occ.addPoint(position[0]+size[1]*np.cos(position[2] * np.pi / 180), position[1]+size[1]*np.sin(position[2] * np.pi / 180), 0, tag=5) # Control mesh size on boundaries
        gmsh.model.occ.addEllipse(position[0], position[1], 0, size[0], size[1], tag=5, zAxis=zAxis, xAxis=xAxis)
        gmsh.model.occ.addCurveLoop([5], tag=2)
        gmsh.model.occ.addPlaneSurface([1, 2], tag=1)
        gmsh.model.occ.mesh.setSize(gmsh.model.occ.getEntities(0), mesh_size)
    
    gmsh.model.occ.synchronize()

    gmsh.model.addPhysicalGroup(1, [2], 1, "inflow")
    gmsh.model.addPhysicalGroup(1, [3], 2, "outflow")
    gmsh.model.addPhysicalGroup(1, [1, 4], 3, "walls")
    gmsh.model.addPhysicalGroup(1, [5], 4, type)
    gmsh.model.addPhysicalGroup(2, [1], 5, "flow_domain")

    
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize('Netgen', niter=5)

    gmsh.write(os.path.join(directory, filename + '.msh'))
    gmsh.finalize()

    msh2xdmf(filename + '.msh', 2, directory)

    np.savez(os.path.join(directory, filename + '.npz'), domain_size=domain_size, size=size, position=position)
    

def createChannelFlowMesh(domain_size, mesh_size = 0.2, filename = '', directory = ''):
    '''
    Generate mesh
    Input:
        domain_size: Absolute size in scale (relative to mesh size). For a nozzle the domain size passes the width of the nozzle as the first element, pass the ratio of expansion as the second element, length of the nozzle as the third element.
            the length of inlet channel and outlet channel are set by default to 5 times the width of the channel for flow stability.
        mesh_size: Mesh size
        filename: Name of the mesh file
        directory: Directory where the mesh file is saved
    Output:
        None
    '''
    gmsh.initialize()
    gmsh.clear()
    gmsh.model.add("t3")

    # Initialize simulation domain
    gmsh.model.occ.addPoint(0, 0, 0, tag=1)
    gmsh.model.occ.addPoint(7 * domain_size[0] + domain_size[2], 0, 0, tag=2)
    gmsh.model.occ.addPoint(0, domain_size[0], 0, tag=3)
    gmsh.model.occ.addPoint(7 * domain_size[0] + domain_size[2], domain_size[0], 0, tag=4)
    gmsh.model.occ.addPoint(2 * domain_size[0], domain_size[0], 0, tag=5)
    gmsh.model.occ.addPoint(2 * domain_size[0] + domain_size[2], domain_size[0] / domain_size[1], 0, tag=6)
    gmsh.model.occ.addPoint(2 * domain_size[0] + domain_size[2], domain_size[0], 0, tag=7)
 
    gmsh.model.occ.addLine(1, 2, 1)
    gmsh.model.occ.addLine(1, 3, 2)
    gmsh.model.occ.addLine(2, 4, 3)
    gmsh.model.occ.addLine(3, 5, 4)
    gmsh.model.occ.addLine(5, 6, 5)
    gmsh.model.occ.addLine(6, 7, 6)
    gmsh.model.occ.addLine(7, 4, 7)
    gmsh.model.occ.addCurveLoop([1, 2, 3, 4, 5, 6, 7], tag=1)

    gmsh.model.occ.addPlaneSurface([1], tag=1)
    gmsh.model.occ.mesh.setSize(gmsh.model.occ.getEntities(0), mesh_size)

    gmsh.model.occ.synchronize()

    gmsh.model.addPhysicalGroup(1, [2], 1, "inflow")
    gmsh.model.addPhysicalGroup(1, [3], 2, "outflow")
    gmsh.model.addPhysicalGroup(1, [1, 4, 5, 6, 7], 3, "walls")
    gmsh.model.addPhysicalGroup(2, [1], 4, "flow_domain")

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize('Netgen', niter=5)

    gmsh.write(os.path.join(directory, filename + '.msh'))
    gmsh.finalize()

    msh2xdmf(filename + '.msh', 2, directory)

    np.savez(os.path.join(directory, filename + '.npz'), domain_size=domain_size)

if __name__ == "__main__":
    
    domain_size = np.array([20, 10])

    # generate 1000 samples for cylinders
    for i in range(10):
        x_coordinate = random.uniform(4, 10)
        y_coordinate = random.uniform(2, 8)
        r = random.uniform(1, np.min([x_coordinate-3, domain_size[0]-x_coordinate-1, y_coordinate-1, domain_size[1]-y_coordinate-1]))
        # createFreeFlowMesh(domain_size, 'circle', np.array([x_coordinate, y_coordinate]), r, mesh_size=0.3, filename='circle_{}'.format(i), directory='io_operations/data221120/has')
        # createFreeFlowMesh(domain_size, 'circle', np.array([x_coordinate, y_coordinate]), r, mesh_size=2, filename='circle_{}'.format(i), directory='io_operations/data221120/las')

    # generate 1000 samples for ellipses
    for j in range(10):
        x_coordinate = random.uniform(4, 10)
        y_coordinate = random.uniform(2, 8)
        a = random.uniform(1, np.min([x_coordinate-3, domain_size[0]-x_coordinate-1, y_coordinate-1, domain_size[1]-y_coordinate-1]))
        b = random.uniform(0.5, 0.7 * np.min([x_coordinate-3, domain_size[0]-x_coordinate-1, y_coordinate-1, domain_size[1]-y_coordinate-1]))
        long_axis = np.max([a, b])
        short_axis = np.min([a, b])
        angle = random.uniform(0, 180)
        # createFreeFlowMesh(domain_size, 'ellipse', np.array([x_coordinate, y_coordinate, angle]), np.array([long_axis, short_axis]), mesh_size=0.3, filename='ellipse_{}'.format(j), directory='io_operations/data221129/has')
        # createFreeFlowMesh(domain_size, 'ellipse', np.array([x_coordinate, y_coordinate, angle]), np.array([long_axis, short_axis]), mesh_size=2, filename='ellipse_{}'.format(j), directory='io_operations/data221129/las')

    # generate 1000 samples for channels
    for k in range(10):
        channel_width = 5
        nozzle_length = random.uniform(0.5, 6)
        expansion_ratio = random.uniform(2, 5)

        domain_size = np.array([channel_width, expansion_ratio, nozzle_length])
        createChannelFlowMesh(domain_size, mesh_size=0.1, filename='channel_{}'.format(k), directory='io_operations/data221206/has')
        createChannelFlowMesh(domain_size, mesh_size=0.5, filename='channel_{}'.format(k), directory='io_operations/data221206/las')