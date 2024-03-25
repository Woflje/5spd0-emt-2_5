import numpy as np

import Integrator_bastest as EJHM
import integrator_bastest_Duffy as EJHM_Duffy
import Basic_Integrator as simple_int
import EFIE_functions as EF
from Mesh import Mesh

from E_field import E_field_in_position

from FastNP import fastNorm

from Timer import Timer

from parameters import parameters
from E_field import E_field_essentials
from integrators.dunavant.Dunavant import get_Dunavant_values

from loose import restructure_mesh_rwg_data

from Integration import int_test

# Initialize the constants
k = 2*np.pi/parameters["wavelength"]

# Initialize dunavant values
dunavant_values, dunavant_weight = get_Dunavant_values(parameters["order_dunavant"])

# Intialize the integrators for the given wavelength
dunavant = EJHM.integrator_bastest(k,parameters["order_dunavant"],dunavant_weight)
duffy = EJHM_Duffy.integrator_bastest(dunavant,k,parameters["order_duffy"],parameters["order_dunavant"],dunavant_values,dunavant_weight)

# Initialize the code which creates the mesh of the input object and sorts all edges
mesh = Mesh(parameters["input_mesh_path"], True)
mesh.prepare_plot_data()  # Always run .plot() because the edge sorting is done here
# mesh.plot_objects()

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

n = 0
while n<N:
    E[n] = int_test(Efield_in,vertices[rwgs[n,0:3]],areas[rwgs[n,4]],dunavant_positions[rwgs[n,4]],dunavant_weight)\
          - int_test(Efield_in,vertices[rwgs[n,[0,1,3]]],areas[rwgs[n,5]],dunavant_positions[rwgs[n,5]],dunavant_weight) #integral E(r)*t(r) dr over surface triangle
    n=n+1

# Plot the incident electric field over the inner edges
mesh.plot_current(E,e_vertex,parameters["E_field_in"]["direction"],parameters["E_field_in"]["polarization"])

# Create system Matrix
A = (N,N)
A = np.zeros(A,dtype=np.complex128)

n=0
i=0

while n<N:  # Loop through all basis functions
    while i<N:  # Loop through all test functions
        if i == n or i > n:  # This is an approximation but saves a lot of computation time (see comment below)*
            for t1 in range(2,4):
                for t2 in range(2,4):
                    factor = -1 if t1 != t2 else 1
                    if EF.check_sing(vertices[rwgs[n,0]],vertices[rwgs[n,1]],vertices[rwgs[n,t1]],vertices[rwgs[i,0]],vertices[rwgs[i,1]],vertices[rwgs[i,t2]]):
                        A[n,i] = A[n,i]+factor*4*np.pi*duffy.int_bastest(vertices[rwgs[n,[0,1,t1]]],vertices[rwgs[i,[0,1,t2]]])
                    else:
                        A[n,i] = A[n,i]+factor*4*np.pi*dunavant.int_bastest(
                            vertices[rwgs[n,[0,1,t1]]],
                            vertices[rwgs[i,[0,1,t2]]],
                            areas[rwgs[n,t1+2]],
                            areas[rwgs[i,t2+2]],
                            dunavant_positions[rwgs[n,t1+2]],
                            dunavant_positions[rwgs[i,t2+2]]
                            )
        i=i+1
    n=n+1
    i=0

A = A+np.transpose(A)-np.diag(np.diag(A)) #this is an approximation but saves a lot of computation time (see comment below)*

# * Only the upper triangle of the system matrix is calculated and is than mirrored. This is an approximation and should be removed for more precision, however it saves a lot of computation time


# Solve the system to find the surface current density
J = np.dot(np.linalg.inv(A),E)

# Plot the currents over the inner edges
mesh.plot_current(J, e_vertex, parameters["E_field_in"]["direction"], parameters["E_field_in"]["polarization"])
#%%
# Start of the post-processing
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
EF.plot_phi(parameters["E_farfield"]["direction"][0],E_farfield)
EF.plot_theta(parameters["E_farfield"]["direction"][1],E_farfield)

# Finalization of the code and delete of all integrator instances
del dunavant
del duffy