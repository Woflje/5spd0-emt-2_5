import numpy as np

from integrators.dunavant.integrator_base_test import Integrator_base_test_dunavant
from integrators.integrator_base_test_duffy import Integrator_base_test_duffy
from integrators.dunavant.dunavant import get_dunavant_values
from integrators.integration import int_test
from plot import plot_phi, plot_theta
from mesh.Mesh import Mesh, check_triangle_pair_singularity, restructure_mesh_rwg_data
from E_field import E_field_in_position, E_field_essentials
from Timer import Timer

from parameters import parameters

# Initialize the constants
k = 2*np.pi/parameters["wavelength"]

# Initialize dunavant values
dunavant_values, dunavant_weight = get_dunavant_values(parameters["order_dunavant"])

# Intialize the integrators for the given wavelength
dunavant = Integrator_base_test_dunavant(k,dunavant_weight)
duffy = Integrator_base_test_duffy(k,parameters["order_Gauss_duffy"],dunavant_values,dunavant_weight)

# Initialize the code which creates the mesh of the input object and sorts all edges
mesh = Mesh(parameters["input_mesh_path"])
mesh.prepare_plot_data()  # Always run .plot() because the edge sorting is done here
if parameters["plots"]:
    mesh.plot_objects()

# Load all edges in r_vect with a xyz for each of the 3 vertices of a Triangle
#  r_vect has N elements for each inner edge with 4 vertices: 0 and 1 for the inner edge and 2 and 3 for the two other vertices that make the 2 triangles

[length, e_vertex, other_vertex, area] = mesh.getmatrix()
N = len(e_vertex)

r_vect = np.array([e_vertex[:, :3], e_vertex[:, 3:], other_vertex[:, :3], other_vertex[:, 3:]]).transpose(1, 0, 2)

vertices, rwgs, areas, dunavant_positions = restructure_mesh_rwg_data(r_vect, dunavant_values)

[E_in_k_vector, E_in_polarization] = E_field_essentials(k, parameters["E_field_in"]["direction"], parameters["E_field_in"]["polarization"])

E = np.zeros([N,1],dtype=np.complex128)

def Efield_in(pos):
    return E_field_in_position(pos, E_in_k_vector, E_in_polarization, parameters["E_field_in"]["amplitude"])

#  Integrate the incident electric field over all test functions

for n in range(N):
    E[n] = int_test(Efield_in,vertices[rwgs[n,0:3]],areas[rwgs[n,4]],dunavant_positions[rwgs[n,4]],dunavant_weight)\
          - int_test(Efield_in,vertices[rwgs[n,[0,1,3]]],areas[rwgs[n,5]],dunavant_positions[rwgs[n,5]],dunavant_weight) #integral E(r)*t(r) dr over surface triangle

# Plot the incident electric field over the inner edges
if parameters["plots"]:
    mesh.plot_current(E,e_vertex,parameters["E_field_in"]["direction"],parameters["E_field_in"]["polarization"])

# Create system Matrix
A = (N,N)
A = np.zeros(A,dtype=np.complex128)

singularities_map = check_triangle_pair_singularity(rwgs)

print('Starting integration')
with Timer('Integration'):
    for n in range(N):
        for i in range(n,N):
            for t1 in range(2,4):
                for t2 in range(2,4):
                    factor = 1 if t1 == t2 else -1
                    if singularities_map[n,i,2*t1+t2-6]:
                        A[n,i] = A[n,i]+factor*4*np.pi*duffy.integrate_base_test(
                            vertices[rwgs[n,[0,1,t1]]],
                            vertices[rwgs[i,[0,1,t2]]],
                            areas[rwgs[n,t1+2]],
                            areas[rwgs[i,t2+2]],
                            dunavant_positions[rwgs[n,t1+2]]
                            )
                    else:
                        A[n,i] = A[n,i]+factor*4*np.pi*dunavant.integrate_base_test(
                            vertices[rwgs[n,[0,1,t1]]],
                            vertices[rwgs[i,[0,1,t2]]],
                            areas[rwgs[n,t1+2]],
                            areas[rwgs[i,t2+2]],
                            dunavant_positions[rwgs[n,t1+2]],
                            dunavant_positions[rwgs[i,t2+2]]
                            )

A = A+np.transpose(A)-np.diag(np.diag(A)) #this is an approximation but saves a lot of computation time (see comment below)*

# * Only the upper triangle of the system matrix is calculated and is than mirrored. This is an approximation and should be removed for more precision, however it saves a lot of computation time

# Solve the system to find the surface current density
J = np.dot(np.linalg.inv(A),E)

# Plot the currents over the inner edges
mesh.plot_current(J, e_vertex, parameters["E_field_in"]["direction"], parameters["E_field_in"]["polarization"])
#%%
# Start of the post-processing
with Timer('Farfield processing'):
    n = 0
    E_ff = np.zeros([N,1],dtype=np.complex128)
    E_farfield = np.zeros((1,np.size(parameters["E_farfield"]["direction"][0]), np.size(parameters["E_farfield"]["direction"][1]), np.size(parameters["E_farfield"]["polarization"])),dtype=np.complex128)

    # For reach range of theta, phi or polarization calculate the farfield Efield
    for a in range(np.size(parameters["E_farfield"]["direction"][0])):
        for b in range(np.size(parameters["E_farfield"]["direction"][1])):
            for c in range(np.size(parameters["E_farfield"]["polarization"])):
                # Calculate the farfield contribution of all edges
                [E_farfield_k_vector, E_farfield_polarization] = E_field_essentials(k, [parameters["E_farfield"]["direction"][0][a],parameters["E_farfield"]["direction"][1][b]], parameters["E_farfield"]["polarization"][c])
                def Efield_ff(pos):
                    return E_field_in_position(pos, E_farfield_k_vector, E_farfield_polarization, parameters["E_farfield"]["amplitude"])
                n = 0
                E_ff = np.zeros([N,1],dtype=np.complex128)
                while n<N:
                    E_ff[n] = J[n]*(int_test(Efield_ff,vertices[rwgs[n,0:3]],areas[rwgs[n,4]],dunavant_positions[rwgs[n,4]],dunavant_weight) - int_test(Efield_ff,vertices[rwgs[n,[0,1,3]]],areas[rwgs[n,5]],dunavant_positions[rwgs[n,5]],dunavant_weight))  #integral E(r)*t(r) dr over surface triangle
                    n=n+1

                # Sum all farfield contributions for a single location
                E_farfield[0,a,b,c] = np.sum(E_ff)

# Plot the farfield results for a constant phi and a constant theta (the input parameters will decide which one is plotted)
plot_phi(parameters["E_farfield"]["direction"][0],E_farfield)
plot_theta(parameters["E_farfield"]["direction"][1],E_farfield)

# Finalization of the code and delete of all integrator instances
del dunavant
del duffy