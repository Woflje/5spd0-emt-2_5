###############

#verify the integrals on a predefined geometry

##############
# Import numpy, the integrator codes, extra functions and the mesh code
import numpy as np
import Integrator_bastest as EJHM
import integrator_bastest_Duffy as EJHM_Duffy
import Basic_Integrator as simple_int
import EFIE_functions as EF
from Mesh import Mesh

############ SET INPUT PARAMETERS###################
wavelength = 1/5
input_mesh = "examples/tetrahedron.dat"
input_EFIELD_AMPLITUDE = 1
input_EFIELD_DIRECTION = [np.pi, np.pi/2]  #angle: phi, theta
input_EFIELD_POLARIZATION = -np.pi/2  #angle of polarization relative to phi and theta

#input_FARFIELD_DIRECTION = np.array([np.linspace(0.001,2*np.pi,180), np.array([np.pi/4.])]) #angle: phi, theta
input_FARFIELD_DIRECTION = np.array([np.array([np.pi/4.]), np.linspace(0.001,2*np.pi,180)]) #angle: phi, theta
input_FARFIELD_POLARIZATION = np.array([np.pi/4.]) #angle of polarization relative to phi and theta


# Initialize the constants
k = 2*np.pi/wavelength


# Intialize the integrators for the given wavelength
dunavant = EJHM.integrator_bastest(k,5)
simple = simple_int.Basic_Integrator(5)
duffy = EJHM_Duffy.integrator_bastest(k,5,5)




# Initialize the code which creates the mesh of the input object and sorts all edges
#mesh = Mesh(input_mesh)
#mesh.plot()  # Always run .plot() because the edge sorting is done here
#[length, e_vertice, other_vertice, area] = mesh.getmatrix()
N = 2

# Load all edges in r_vect with a xyz for each of the 3 vertices of a Triangle
#  r_vect has N elements for each inner edge with 4 vertices: 0 and 1 for the inner edge and 2 and 3 for the two other vertices that make the 2 triangles
r_vect = np.empty([N,4,3])  # Position vectors of the basis and test functions
#n = 0
#while n<N:
#    r_vect[n] = np.array([np.array([e_vertice[n,0], e_vertice[n,1], e_vertice[n,2]]), np.array([e_vertice[n,3], e_vertice[n,4], e_vertice[n,5]]), np.array([other_vertice[n,0], other_vertice[n,1], other_vertice[n,2]]), np.array([other_vertice[n,3], other_vertice[n,4], other_vertice[n,5]])])
#    n = n+1  #update within array of vertices
r_vect[0]=np.array([np.array([0.0038,0.0026,-0.0008]),np.array([-0.001,0.0022,0.002]),np.array([0.0024,0.0013,-0.0016]),np.array([0.0018,0.005,0.003])])
r_vect[1]=np.array([np.array([0.014,0.0024,-0.004]),np.array([0.008,0.002,0]),np.array([0.012,-0.003,0.004]),np.array([0.021,0.009,-0.009])])
# Create and integrate the incomming Efield
E = np.zeros([N,1],dtype=np.complex128)
def Efield_in(pos):
    return EF.E_in(input_EFIELD_AMPLITUDE, input_EFIELD_DIRECTION, pos, input_EFIELD_POLARIZATION, wavelength)

#  Integrate the incident electric field over all test functions
n = 0
while n<N:
    E[n] = simple.int_test(Efield_in,np.array([r_vect[n,0], r_vect[n,1], r_vect[n,2]])) - simple.int_test(Efield_in,np.array([r_vect[n,0], r_vect[n,1], r_vect[n,3]])) #integral E(r)*t(r) dr over surface triangle
    n=n+1

# Plot the incident electric field over the inner edges
#mesh.plot_current(E,e_vertice,input_EFIELD_DIRECTION,input_EFIELD_POLARIZATION)

# Create system Matrix
A = (N,N)
A = np.zeros(A,dtype=np.complex128)

n=0
i=0
while n<N:  # Loop through all basis functions
    while i<N:  # Loop through all test functions
        if i == n or i > n:  # This is an approximation but saves a lot of computation time (see comment below)*
            # First check for singularity to decide which integrator to use and than integrate each triangle of each RWG over each other
            print("\n")
            if EF.check_sing(r_vect[n, 0], r_vect[n, 1], r_vect[n, 2], r_vect[i, 0], r_vect[i, 1], r_vect[i, 2]):
                int_temp = duffy.int_bastest(np.array([r_vect[n,0], r_vect[n,1], r_vect[n,2]]), np.array([r_vect[i,0], r_vect[i,1], r_vect[i,2]]))
                A[n,i] = A[n,i] + int_temp
                print("duffy",n,i,int_temp,"\n sum:",A[n,i])
            else:
                int_temp = dunavant.int_bastest(np.array([r_vect[n, 0], r_vect[n, 1], r_vect[n, 2]]), np.array([r_vect[i, 0], r_vect[i, 1], r_vect[i, 2]]))
                A[n, i] = A[n, i] + int_temp
                print("dunavant",n,i,int_temp,"\n sum:",A[n,i])
            if EF.check_sing(r_vect[n,0], r_vect[n,1], r_vect[n,3], r_vect[i,0], r_vect[i,1], r_vect[i,2]):
                int_temp = duffy.int_bastest(np.array([r_vect[n,0], r_vect[n,1], r_vect[n,3]]), np.array([r_vect[i,0], r_vect[i,1], r_vect[i,2]]))
                A[n,i] = A[n,i] - int_temp
                print("duffy",n,i,int_temp,"\n sum:",A[n,i])
            else:
                int_temp = dunavant.int_bastest(np.array([r_vect[n,0], r_vect[n,1], r_vect[n,3]]), np.array([r_vect[i,0], r_vect[i,1], r_vect[i,2]]))
                A[n,i] = A[n,i] - int_temp
                print("dunavant",n,i,int_temp,"\n sum:",A[n,i])
            if EF.check_sing(r_vect[n,0], r_vect[n,1], r_vect[n,2], r_vect[i,0], r_vect[i,1], r_vect[i,3]):
                int_temp = duffy.int_bastest(np.array([r_vect[n, 0], r_vect[n, 1], r_vect[n, 2]]), np.array([r_vect[i, 0], r_vect[i, 1], r_vect[i, 3]]))
                A[n,i] = A[n,i] - int_temp
                print("duffy",n,i,int_temp,"\n sum:",A[n,i])
            else:
                int_temp = dunavant.int_bastest(np.array([r_vect[n, 0], r_vect[n, 1], r_vect[n, 2]]), np.array([r_vect[i, 0], r_vect[i, 1], r_vect[i, 3]]))
                A[n, i] = A[n, i] - int_temp
                print("dunavant",n,i,int_temp,"\n sum:",A[n,i])
            if EF.check_sing(r_vect[n,0], r_vect[n,1], r_vect[n,3], r_vect[i,0], r_vect[i,1], r_vect[i,3]):
                int_temp = duffy.int_bastest(np.array([r_vect[n, 0], r_vect[n, 1], r_vect[n, 3]]), np.array([r_vect[i, 0], r_vect[i, 1], r_vect[i, 3]]))
                A[n, i] = A[n, i] + int_temp
                print("duffy",n,i,int_temp,"\n result:",A[n,i])
            else:
                int_temp = dunavant.int_bastest(np.array([r_vect[n, 0], r_vect[n, 1], r_vect[n, 3]]), np.array([r_vect[i, 0], r_vect[i, 1], r_vect[i, 3]]))
                A[n, i] = A[n, i] + int_temp
                print("dunavant",n,i,int_temp,"\n result:",A[n,i])
        i=i+1
    n=n+1
    i=0

#Check het optellen
#Check de angles

A = (A+np.transpose(A)-np.diag(np.diag(A))) #this is an approximation but saves a lot of computation time (see comment below)*
print("\n",A)
# * Only the upper triangle of the system matrix is calculated and is than mirrored. This is an approximation and should be removed for more precision, however it saves a lot of computation time


# Solve the system to find the surface current density
J = np.dot(np.linalg.inv(A),E)

# Plot the currents over the inner edges
#mesh.plot_current(J, e_vertice, input_EFIELD_DIRECTION, input_EFIELD_POLARIZATION)

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
#EF.plot_phi(input_FARFIELD_DIRECTION[0],E_farfield)
#EF.plot_theta(input_FARFIELD_DIRECTION[1],E_farfield)

# Finalization of the code and delete of all integrator instances
del dunavant
del duffy
del simple
