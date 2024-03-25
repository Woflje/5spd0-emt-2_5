import numpy as np
from Mesh import Mesh

wavelength = 100
input_mesh = "examples/tetrahedron.dat"

# E-Field Amplitude, Direction [phi,theta], Polarization angle relative to phi and theta
input_EFIELD_AMPLITUDE = 1
input_EFIELD_DIRECTION = [np.pi/4,np.pi/4]
input_EFIELD_POLARIZATION = np.pi/4

input_FARFIELD_DIRECTION = [np.array([np.pi*2.]), np.linspace(0.001,2*np.pi,180)] #angle: phi, theta
input_FARFIELD_POLARIZATION = np.array([0]) #angle of polarization relative to phi and theta

# Initialize the constants
k = 2*np.pi/wavelength
def fastCross(v1,v2):
    #reshape vectors into shape accepted by np.einsum
    v1=np.reshape(v1,(1,3))
    v2=np.reshape(v2,(1,3))
    #use einstein sum notation to perform cross product
    eijk=eijk = np.zeros((3, 3, 3))
    eijk[0, 1, 2] = eijk[1, 2, 0] = eijk[2, 0, 1] = 1
    eijk[0, 2, 1] = eijk[2, 1, 0] = eijk[1, 0, 2] = -1
    s=np.einsum('iuk,vk->uvi', np.einsum('ijk,uj->iuk', eijk, v1), v2)
    return np.reshape(s,(3)) #reshape result so the existing fucntions accept it

def fastNorm(v):
    return np.sqrt(np.sum(v**2, axis=-1))

def int_triangle(self,Phi,T):
    A_tot = area_triangle(T) #determine area of triangle
    DunPolys = Dunavant_Values.PolyDun() #Initailize function to obtain dunnavant values
    weight, alpha, beta, gamma, ng, rot = DunPolys.Dunavant_Values(self.order) #obtain the dunvanat values of the given order
    #temp = 0;
    
    #Integral over triangle T with a scalar function according to dunavant method in barycentric coordinates.
    r_v=np.dot(np.vstack([alpha,beta,gamma]).transpose(),T)
    phi_t=[Phi(r) for r in r_v]
    temp=sum(np.multiply(phi_t,weight))
    return A_tot*temp

def int_test(E,T):    
    l_edge = fastNorm(np.subtract(T[0],T[1])) #Length of common edge
    A = area_triangle(T)
    scalar = lambda r: np.dot(l_edge/(2*A)*np.subtract(r,T[2]),E(r)) #This function ensures that we only have to integrate over a scalar function as the test/basis function dotted with the vector field results in a scalar function.
    
    Itest = int_triangle(scalar,T) #Integrate the new scalar function using the int_triangle function.
    
    return Itest

def E_in(amp, dir, pos, polarization, wavelength): # Function to calculate the incoming electric field
    #Determine wavenumber
    k = 2 * np.pi / (wavelength)
    #Determine wave vector using input angles phi (dir[0]) and theta (dir[1])
    kx = k*np.sin(dir[0])*np.cos(dir[1])
    ky = k*np.sin(dir[0])*np.sin(dir[1])
    kz = k*np.cos(dir[0])

    #Find values for a vector with a slightly different phi value
    kx_phi = k*np.sin(dir[0]+0.001)*np.cos(dir[1])
    ky_phi = k*np.sin(dir[0]+0.001)*np.sin(dir[1])
    kz_phi = k*np.cos(dir[0]+0.001)
    #Cross product of the wave vector and the changed phi vector to find a perpendicular vector
    Vect1 = fastCross(np.array([kx_phi,ky_phi,kz_phi]), np.array([kx,ky,kz]))

    #Find values for a vector with a slightly different theta value
    kx_th = k*np.sin(dir[0])*np.cos(dir[1]+0.004)
    ky_th = k*np.sin(dir[0])*np.sin(dir[1]+0.004)
    kz_th = k*np.cos(dir[0])
    #Cross product of the wave vector and the changed theta vector to find a perpendicular vector
    Vect2 = fastCross(np.array([kx_th,ky_th,kz_th]), np.array([kx,ky,kz])) 

    #Normalize the perpendicular vectors
    pol_vect1 = np.multiply((1./fastNorm(Vect1)),Vect1)
    #print('k phi')
    #print(pol_vect1)
    pol_vect2 = np.multiply((1./fastNorm(Vect2)),Vect2)
    #print('k theta')
    #print(pol_vect2)
    #Allow user to choose polarization by using input angle between the two normalized polarization vectors
    polarization = np.sin(polarization)*pol_vect2+np.cos(polarization)*pol_vect1
    polarization = np.multiply((1/fastNorm(polarization)),polarization)
    #print('polarization')
    #print(polarization)

    #Check if polarization is perpendicular to the wave vector
    if np.dot([kx, ky, kz], polarization) > 1e-10:
        print("The polarization and wavenumber are not perpendicular!")

    #return the full plane wave
    #return np.array([1,0,0])
    return np.multiply(amp*polarization,np.exp(np.multiply(-1J, np.dot([kx, ky, kz],pos))))



# Initialize the code which creates the mesh of the input object and sorts all edges
mesh = Mesh(input_mesh, True)
mesh.prepare_plot_data()

[length, e_vertex, other_vertex, area] = mesh.getmatrix()
N = len(e_vertex)

r_vect = np.array([e_vertex[:, :3], e_vertex[:, 3:], other_vertex[:, :3], other_vertex[:, 3:]]).transpose(1, 0, 2)

# Create and integrate the incomming Efield
E= np.zeros([N,1],dtype=np.complex128)

def Efield_in(pos):
    return E_in(input_EFIELD_AMPLITUDE, input_EFIELD_DIRECTION, pos, input_EFIELD_POLARIZATION, wavelength)

#  Integrate the incident electric field over all test functions

n = 0
while n<N:
    E[n] = int_test(Efield_in,np.array([r_vect[n,0], r_vect[n,1], r_vect[n,2]])) - int_test(Efield_in,np.array([r_vect[n,0], r_vect[n,1], r_vect[n,3]])) #integral E(r)*t(r) dr over surface triangle
    n=n+1