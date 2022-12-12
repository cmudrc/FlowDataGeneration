from dolfin import *
import numpy as np
from msh2xdmf import import_mesh
import matplotlib.pyplot as plt

predicted_dir = './data/predicted/data_4999.npz'

mesh, mf_boundaries, association_table = import_mesh(prefix="cylinder_las", subdomains=False, directory="io_operations")

gdim = mesh.geometry().dim()
nv = mesh.num_vertices()

def mesh2triang(mesh):
    import matplotlib.tri as tri
    xy = mesh.coordinates()
    return tri.Triangulation(xy[:, 0], xy[:, 1], mesh.cells())

predicted_result = np.load(predicted_dir)
u = predicted_result['u']
p = predicted_result['p']

X = mesh.coordinates()
X = [X[:, i] for i in range(gdim)]
U = [u[i * nv: (i + 1) * nv] for i in range(gdim)]

ax = plt.gca()
ax.set_aspect('equal')

C = U[0]**2
for i in range(1, gdim):
    C += U[i]**2
C = np.sqrt(C)

args = X + U + [C]
ax.quiver(*args)

plt.figure()
ax.tricontourf(mesh2triang(mesh), p, 40)
plt.show()