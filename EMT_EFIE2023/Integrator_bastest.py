import numpy as np
from FastNP import fastNorm

from Integration import int_triangle


class integrator_bastest():
    def __init__(self,k,order,dunavant_weight):
        #Initialize integrator
        self.mu = 4*np.pi*10**(-7) #permeability of free space
        self.c = 299792458 #light speed
        self.ep = 1/(self.mu*self.c**2) #permitivity of free space
        self.omega = self.c*k #radial frequency
        self.k = k # wave vector
        self.order = order #order of the dunavant method
        self.dunavant_weight = dunavant_weight
    
    def int_bastest(self,T1,T2,A1,A2,dunavant_pos1,dunavant_pos2): #initialize the basic integrators
        length1 = fastNorm(np.subtract(T1[1],T1[0])) #length common edge
        length2 = fastNorm(np.subtract(T2[1],T2[0])) #lengthe common edge 
        j = complex(0,1) #define imaginairy unit 
        
        norm = lambda r,rp: fastNorm(np.subtract(r,rp)) #function to calculate distance between observation point and source point
        Scalar_green = lambda r,rp: np.exp(-j*self.k*norm(r,rp))/(4*np.pi*norm(r,rp)) #function for the scalar green function  
        Term1 = lambda r,rp: length1/(2*A1)*(r-T1[2]) #test function
        Term2 = lambda r,rp: length2/(2*A2)*(rp-T2[2]) #basis function
        
        Part1 = lambda r,rp: 1/(j*self.k)*length2/A2*length1/A1*Scalar_green(rp,r) #divergence of test and basis function with scalar green
        Part2 = lambda r,rp: j*self.k*np.dot(Term1(r,rp),Term2(r,rp))*Scalar_green(rp,r) #basis test function and scalar green
     
        Int1 = lambda r: int_triangle(lambda rp: Part1(r,rp),A2,dunavant_pos2,self.dunavant_weight) #Hyper singular integral calculation
        Int2 = lambda r: int_triangle(lambda rp: Part2(r,rp),A2,dunavant_pos2,self.dunavant_weight) #singular integral calculation
        Int_result = int_triangle(Int1,A1,dunavant_pos1,self.dunavant_weight) + int_triangle(Int2,A2,dunavant_pos2,self.dunavant_weight) #Total result is the sum of the singular and hyper singular integrals

        return Int_result