import numpy as np
from scipy.io import savemat

from integrators.dunavant.integrator_base_test import Integrator_base_test_dunavant
from integrators.integrator_base_test_duffy import Integrator_base_test_duffy
from integrators.dunavant.dunavant import get_dunavant_values

from plot import plot_phi, plot_theta
from mesh.Mesh import Mesh, check_triangle_pair_intersections, restructure_mesh_rwg_data
from E_field import E_field
from Timer import Timer, timestamp

from parameters import parameters

timestamp('Start code')

k = 2*np.pi/parameters["wavelength"]

dunavant_values, dunavant_weight = get_dunavant_values(parameters["order_dunavant"])

mesh = Mesh(parameters["input_mesh_path"])
mesh.prepare_plot_data() # This sorts the data

[length, e_vertex, other_vertex, area] = mesh.getmatrix()
N = len(e_vertex) # N is the number of RWGs
r_vect = np.array([e_vertex[:, :3], e_vertex[:, 3:], other_vertex[:, :3], other_vertex[:, 3:]]).transpose(1, 0, 2)
vertices, rwgs, areas, dunavant_positions = restructure_mesh_rwg_data(r_vect, dunavant_values)

E_incident = E_field(vertices, rwgs, areas, dunavant_positions, k, parameters["E_field_in"]["amplitude"], parameters["E_field_in"]["direction"], parameters["E_field_in"]["polarization"], dunavant_weight)

# Plot the incident electric field over the inner edges
if parameters["plots"]:
    mesh.plot_current(E_incident,e_vertex,parameters["E_field_in"]["direction"],parameters["E_field_in"]["polarization"])

# Create system Matrix
A = (N,N)
A = np.zeros(A,dtype=np.complex128)

dunavant = Integrator_base_test_dunavant(k,dunavant_weight)
duffy = Integrator_base_test_duffy(k,parameters["order_Gauss_duffy"],dunavant_weight)

intersections_map = check_triangle_pair_intersections(rwgs)

# Integration for base and test functions
with Timer('Integration',True):
    for n in range(N):
        for i in range(n,N): # Starting at n, to only compute upper triangle
            for t1 in range(2,4):
                for t2 in range(2,4):
                    factor = 1 if t1 == t2 else -1
                    T1 = rwgs[n,[0,1,t1]]
                    T2 = rwgs[i,[0,1,t2]]
                    if intersections_map[n,i,2*t1+t2-6]:
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

# Only upper triangle including diagonal is computed and mirrored for a time saving approximation
A = A+np.transpose(A)-np.diag(np.diag(A))

# Solve the system to find the surface current density
J = np.dot(np.linalg.inv(A),E_incident)

# Plot the current densities over the inner edges
if parameters["plots"]:
    mesh.plot_current(J, e_vertex, parameters["E_field_in"]["direction"], parameters["E_field_in"]["polarization"])

with Timer('Farfield processing',True):
    E_farfield = np.zeros((1,np.size(parameters["E_farfield"]["direction"][0]), np.size(parameters["E_farfield"]["direction"][1]), np.size(parameters["E_farfield"]["polarization"])),dtype=np.complex128)

    # For each range of theta, phi or polarization calculate the farfield
    for a in range(np.size(parameters["E_farfield"]["direction"][0])):
        for b in range(np.size(parameters["E_farfield"]["direction"][1])):
            for c in range(np.size(parameters["E_farfield"]["polarization"])):
                E_farfield_per_angle = J*E_field(
                    vertices, rwgs, areas, dunavant_positions, k, parameters["E_field_in"]["amplitude"],
                    [parameters["E_farfield"]["direction"][0][a],parameters["E_farfield"]["direction"][1][b]],
                    parameters["E_farfield"]["polarization"][c], dunavant_weight
                )
                # Sum all farfield contributions for a single location
                E_farfield[0,a,b,c] = np.sum(E_farfield_per_angle)

# Plot the farfield results for a constant phi and a constant theta
if parameters["plots"]:
    plot_phi(parameters["E_farfield"]["direction"][0],E_farfield)
    plot_theta(parameters["E_farfield"]["direction"][1],E_farfield)

if parameters["export_E_farfield"]:
    mat_data = {
        parameters["export_E_farfield_data_name"]: E_farfield
    }
    savemat(f"{parameters['export_E_farfield_filename']}.mat",mat_data)