import numpy as np

import Integrator_bastest as EJHM
import integrator_bastest_Duffy as EJHM_Duffy
import Basic_Integrator as simple_int
import EFIE_functions as EF
from Mesh import Mesh

from Timer import Timer

from parameters import parameters
from E_field import E_field_essentials, E_field_in_points
from integrators.dunavant.Dunavant import get_Dunavant_values

from loose import calculate_triangle_area, restructure_mesh_rwg_data

# Initialize the constants
k = 2*np.pi/parameters["wavelength"]

# Initialize dunavant values
dunavant_values, dunavant_weight = get_Dunavant_values(parameters["order_dunavant"])

# Intialize the integrators for the given wavelength
dunavant = EJHM.integrator_bastest(k,5)
duffy = EJHM_Duffy.integrator_bastest(k,5,5)

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

[E_in_k_vector, E_in_polarization] = E_field_essentials(k, parameters["input_E_field"]["direction"], parameters["input_E_field"]["polarization"])

# Compute the electric field at all Dunavant points
E_field_in = E_field_in_points(dunavant_positions, E_in_k_vector, E_in_polarization, parameters["input_E_field"]["amplitude"])

# Vectorized selection of triangle indices for the first and second triangles of each RWG
triangle_indices_first = rwgs[:, 4]
triangle_indices_second = rwgs[:, 5]

# Correct reshaping of dunavant_weight for broadcasting
dunavant_weight_reshaped = dunavant_weight.reshape(1, -1, 1)  # Shape becomes (1, P, 1)

# Recalculate E_contribution_first with the corrected weight shape
E_contribution_first = np.sum(E_field_in[triangle_indices_first] * dunavant_weight_reshaped, axis=(1, 2)) * areas[triangle_indices_first]
E_contribution_second = np.sum(E_field_in[triangle_indices_second] * dunavant_weight_reshaped, axis=(1, 2)) * areas[triangle_indices_second]

# Compute the final E_contribution
E_contribution = E_contribution_first - E_contribution_second

# Reshape E_contribution to match the expected output shape (N, 1)
E_contribution = E_contribution.reshape(-1, 1)
print(1)

def integrate_test(self, integrand, triangle_points, triangle_area):
    l_edge = np.linalg.norm(triangle_points[0]-triangle_points[1]) #Length of common edge
    scalar = lambda r: np.dot(l_edge/(2*triangle_area)*(r-triangle_points[2]),integrand(r)) #This function ensures that we only have to integrate over a scalar function as the test/basis function dotted with the vector field results in a scalar function.
    Itest = self.integrate_known_triangle(scalar,triangle_area) #Integrate the new scalar function using the int_triangle function.
    return Itest

E= np.zeros([N,1],dtype=np.complex128)

def Efield_(pos):
    return EF.E_in(input_EFIELD_AMPLITUDE, input_EFIELD_DIRECTION, pos, input_EFIELD_POLARIZATION, wavelength)


n = 0
for n in range (0,N):
    E[n] = integrate_test(E_field_, np.array(vertices[rwgs[n,0:3]])) - integrate_test(E_field_, np.array([vertices[rwgs[n,0]], vertices[rwgs[n,1]], vertices[rwgs[n,3]]]))

print(1)

# Plot the incident electric field over the inner edges
mesh.plot_current(E,e_vertex,input_EFIELD_DIRECTION,input_EFIELD_POLARIZATION)

approach=0

if approach==0:
    # Create system Matrix
    A = np.zeros((N,N),dtype=np.complex128)
    n=0
    i=0
    while n<N:  # Loop through all basis functions
        while i<N:  # Loop through all test functions
            if i == n or i > n:  # This is an approximation but saves a lot of computation time (see comment below)*
                # First check for singularity to decide which integrator to use and than integrate each triangle of each RWG over each other
                if EF.check_sing(r_vect[n, 0], r_vect[n, 1], r_vect[n, 2], r_vect[i, 0], r_vect[i, 1], r_vect[i, 2]):
                    A[n,i] = A[n,i]+ 4*np.pi*duffy.int_bastest(np.array([r_vect[n,0], r_vect[n,1], r_vect[n,2]]), np.array([r_vect[i,0], r_vect[i,1], r_vect[i,2]]))
                else:
                    A[n, i] = A[n, i] + 4*np.pi*dunavant.int_bastest(np.array([r_vect[n, 0], r_vect[n, 1], r_vect[n, 2]]), np.array([r_vect[i, 0], r_vect[i, 1], r_vect[i, 2]]))
                if EF.check_sing(r_vect[n,0], r_vect[n,1], r_vect[n,2], r_vect[i,0], r_vect[i,1], r_vect[i,3]):
                    A[n,i] = A[n,i] - 4*np.pi*duffy.int_bastest(np.array([r_vect[n,0], r_vect[n,1], r_vect[n,2]]), np.array([r_vect[i,0], r_vect[i,1], r_vect[i,3]]))
                else:
                    A[n,i] = A[n,i] - 4*np.pi*dunavant.int_bastest(np.array([r_vect[n,0], r_vect[n,1], r_vect[n,2]]), np.array([r_vect[i,0], r_vect[i,1], r_vect[i,3]]))
                if EF.check_sing(r_vect[n,0], r_vect[n,1], r_vect[n,3], r_vect[i,0], r_vect[i,1], r_vect[i,2]):
                    A[n,i] = A[n,i] - 4*np.pi*duffy.int_bastest(np.array([r_vect[n, 0], r_vect[n, 1], r_vect[n, 3]]), np.array([r_vect[i, 0], r_vect[i, 1], r_vect[i, 2]]))
                else:
                    A[n, i] = A[n, i] - 4*np.pi*dunavant.int_bastest(np.array([r_vect[n, 0], r_vect[n, 1], r_vect[n, 3]]), np.array([r_vect[i, 0], r_vect[i, 1], r_vect[i, 2]]))
                if EF.check_sing(r_vect[n,0], r_vect[n,1], r_vect[n,3], r_vect[i,0], r_vect[i,1], r_vect[i,3]):
                    A[n, i] = A[n, i] + 4*np.pi*duffy.int_bastest(np.array([r_vect[n, 0], r_vect[n, 1], r_vect[n, 3]]), np.array([r_vect[i, 0], r_vect[i, 1], r_vect[i, 3]]))
                else:
                    A[n, i] = A[n, i] + 4*np.pi*dunavant.int_bastest(np.array([r_vect[n, 0], r_vect[n, 1], r_vect[n, 3]]), np.array([r_vect[i, 0], r_vect[i, 1], r_vect[i, 3]]))
            i=i+1
        n=n+1
        i=0
else:
    A = np.zeros((N,N),dtype=np.complex128)
    is_singular = np.zeros((N,N),dtype=bool)
    for n in range(N):
        for i in range(n, N):
            is_singular[n,i] = EF.check_sing_vectorized(r_vect[n, 0:3], r_vect[i, 0:3])
    for n in range(N):
        for i in range(n, N):
            if is_singular[n, i]:
                A[n, i] += 4 * np.pi * duffy.int_bastest_vectorized(r_vect[n, 0:3], r_vect[i, 0:3])
            else:
                A[n, i] += 4 * np.pi * dunavant.int_bastest_vectorized(r_vect[n, 0:3], r_vect[i, 0:3])
            A[i, n] = A[n, i]  # Mirror result to leverage symmetry

A = A+np.transpose(A)-np.diag(np.diag(A)) #this is an approximation but saves a lot of computation time (see comment below)*

# * Only the upper triangle of the system matrix is calculated and is than mirrored. This is an approximation and should be removed for more precision, however it saves a lot of computation time


# Solve the system to find the surface current density
J = np.dot(np.linalg.inv(A),E)

# Plot the currents over the inner edges
mesh.plot_current(J, e_vertex, input_EFIELD_DIRECTION, input_EFIELD_POLARIZATION)
#%%
# Start of the post-processing
n = 0
E_ff = np.zeros([N,1],dtype=np.complex128)
E_farfield = np.zeros((1,np.size(input_FARFIELD_DIRECTION[0]), np.size(input_FARFIELD_DIRECTION[1]), np.size(input_FARFIELD_POLARIZATION)),dtype=np.complex128)

# For reach range of theta, phi or polarization calculate the farfield Efield
for a in range(np.size(input_FARFIELD_DIRECTION[0])):
    for b in range(np.size(input_FARFIELD_DIRECTION[1])):
        for c in range(np.size(input_FARFIELD_POLARIZATION)):
            # Calculate the farfield contribution of all edges
            def Efield_ff(pos):
                return EF.E_in(input_EFIELD_AMPLITUDE, [input_FARFIELD_DIRECTION[0][a], input_FARFIELD_DIRECTION[1][b]], pos, input_FARFIELD_POLARIZATION[c], wavelength)
            n = 0
            E_ff = np.zeros([N,1],dtype=np.complex128)
            while n<N:
                E_ff[n] = J[n]*(simple.int_test(Efield_ff,np.array([r_vect[n,0], r_vect[n,1], r_vect[n,2]])) - simple.int_test(Efield_ff,np.array([r_vect[n,0], r_vect[n,1], r_vect[n,3]])))  #integral E(r)*t(r) dr over surface triangle
                n=n+1

            # Sum all farfield contributions for a single location
            E_farfield[0,a,b,c] = np.sum(E_ff)


# Plot the farfield results for a constant phi and a constant theta (the input parameters will decide which one is plotted)
EF.plot_phi(input_FARFIELD_DIRECTION[0],E_farfield)
EF.plot_theta(input_FARFIELD_DIRECTION[1],E_farfield)

# Finalization of the code and delete of all integrator instances
del dunavant
del duffy
del simple