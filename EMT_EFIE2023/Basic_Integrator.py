import numpy as np
from FastNP import fastNorm
from integrators.dunavant.Dunavant import get_Dunavant_values

class Basic_Integrator():
    def __init__(self,order):
        self.order = order
        
    #function that calculates the area of a triangle 
    def area_triangle(self,T):
        #lengths of all sides of the given triangle T
        l1 = fastNorm(T[:, 1] - T[:, 0])
        l2 = fastNorm(T[:, 0] - T[:, 2])
        l3 = fastNorm(T[:, 2] - T[:, 1])
        A = 1/4*np.sqrt(4*l1**2*l2**2-(l1**2+l2**2-l3**2)**2)
        return A             
     
    #Function that integrates a scalar green function over a triangle 
    def int_triangle(self,Phi,T,A=None):
        if A is None:
            A = self.area_triangle(T) #determine area of triangle
        weight, alpha, beta, gamma, ng = get_Dunavant_values(order=self.order) #obtain the dunvanat values of the given order
        
        #Integral over triangle T with a scalar function according to dunavant method in barycentric coordinates.
        r_v=np.dot(np.vstack([alpha,beta,gamma]).T,T)
        phi_t = Phi(r_v)
        # phi_t=[Phi(r) for r in r_v]
        return A*sum(np.multiply(phi_t,weight))
    
    #Function that calculates the integral over a triangle for a vector field
    def int_test(self,E,T):    
        l_edge = fastNorm(T[0]-T[1]) #Length of common edge
        A = self.area_triangle(T)
        scalar = lambda r: np.dot(l_edge/(2*A)*(r-T[2]),E(r)) #This function ensures that we only have to integrate over a scalar function as the test/basis function dotted with the vector field results in a scalar function.
        Itest = self.int_triangle(scalar,T,A=A) #Integrate the new scalar function using the int_triangle function.
        return Itest