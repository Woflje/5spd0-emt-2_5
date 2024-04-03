import numpy as np

from integrators.dunavant.integrator_base_test import Integrator_base_test_dunavant
from integrators.integrator_base_test_duffy import Integrator_base_test_duffy
from integrators.dunavant.dunavant import get_dunavant_values
from plot import plot_phi, plot_theta
from mesh.Mesh import Mesh, check_triangle_pair_singularity, restructure_mesh_rwg_data
from E_field import E_field
from Timer import Timer, timestamp

from scipy.io import savemat

from parameters import parameters

timestamp('Start code')

# Initialize the constants
k = 2*np.pi/parameters["wavelength"]

# Initialize dunavant values
dunavant_values, dunavant_weight = get_dunavant_values(parameters["order_dunavant"])

# Intialize the integrators for the given wavelength
dunavant = Integrator_base_test_dunavant(k,dunavant_weight)
duffy = Integrator_base_test_duffy(k,parameters["order_Gauss_duffy"],dunavant_weight)

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

E = E_field(vertices, rwgs, areas, dunavant_positions, k, parameters["E_field_in"]["amplitude"], parameters["E_field_in"]["direction"], parameters["E_field_in"]["polarization"], dunavant_weight)
    

# Plot the incident electric field over the inner edges
if parameters["plots"]:
    mesh.plot_current(E,e_vertex,parameters["E_field_in"]["direction"],parameters["E_field_in"]["polarization"])
# Create system Matrix
A = (N,N)
A = np.zeros(A,dtype=np.complex128)

singularities_map = check_triangle_pair_singularity(rwgs)

with Timer('Integration',True):
    for n in range(N):
        for i in range(n,N):
            for t1 in range(2,4):
                for t2 in range(2,4):
                    factor = 1 if t1 == t2 else -1
                    T1 = rwgs[n,[0,1,t1]]
                    T2 = rwgs[i,[0,1,t2]]
                    if singularities_map[n,i,2*t1+t2-6]:
                        A[n,i] += factor*4*np.pi*duffy.integrate_base_test(
                            vertices[T1],
                            vertices[T2],
                            areas[rwgs[n,t1+2]],
                            areas[rwgs[i,t2+2]],
                            dunavant_positions[n,t1-2]
                            )
                    else:
                        A[n,i] += factor*4*np.pi*dunavant.integrate_base_test(
                            vertices[T1],
                            vertices[T2],
                            areas[rwgs[n,t1+2]],
                            areas[rwgs[i,t2+2]],
                            dunavant_positions[n,t1-2],
                            dunavant_positions[i,t2-2]
                            )

del dunavant
del duffy

A = A+np.transpose(A)-np.diag(np.diag(A)) #this is an approximation but saves a lot of computation time (see comment below)

# Only the upper triangle of the system matrix is calculated and is than mirrored. This is an approximation and should be removed for more precision, however it saves a lot of computation time

# Solve the system to find the surface current density
J = np.dot(np.linalg.inv(A),E)

# Plot the currents over the inner edges
if parameters["plots"]:
    mesh.plot_current(J, e_vertex, parameters["E_field_in"]["direction"], parameters["E_field_in"]["polarization"])

# Start of the post-processing

with Timer('Farfield processing',True):
    n = 0
    E_farfield = np.zeros((1,np.size(parameters["E_farfield"]["direction"][0]), np.size(parameters["E_farfield"]["direction"][1]), np.size(parameters["E_farfield"]["polarization"])),dtype=np.complex128)

    # For reach range of theta, phi or polarization calculate the farfield Efield
    for a in range(np.size(parameters["E_farfield"]["direction"][0])):
        for b in range(np.size(parameters["E_farfield"]["direction"][1])):
            for c in range(np.size(parameters["E_farfield"]["polarization"])):
                E_ff = J*E_field(vertices, rwgs, areas, dunavant_positions, k, parameters["E_field_in"]["amplitude"], [parameters["E_farfield"]["direction"][0][a],parameters["E_farfield"]["direction"][1][b]], parameters["E_farfield"]["polarization"][c], dunavant_weight)

                # Sum all farfield contributions for a single location
                E_farfield[0,a,b,c] = np.sum(E_ff)

# Plot the farfield results for a constant phi and a constant theta (the input parameters will decide which one is plotted)
if parameters["plots"]:
    plot_phi(parameters["E_farfield"]["direction"][0],E_farfield)
    plot_theta(parameters["E_farfield"]["direction"][1],E_farfield)

# Finalization of the code and delete of all integrator instances

E_farfield_output_name = f"E_farfield_matrix_{parameters['input_mesh_path'].split('.')[0].split('/')[-1]}"

mymat = {
    E_farfield_output_name: E_farfield
}
savemat(f"{E_farfield_output_name}.mat",mymat)