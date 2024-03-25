import numpy as np
import scipy
import math
import Basic_Integrator_old as Basic_Integrator
import FastNP as fnp


class integrator_bastest():
    def __init__(self,k,order_duffy,order_dunavant):
        print("Initialize integration Duffy method")
        self.k = k # wave vector
        self.order_duffy = order_duffy #oder of the duffy method
        self.order_dunavant = order_dunavant #order of the dunvant method
        
        #Function that computes the integral over one subdomain. If called three times the whole triangle is integrated using duffy
    def EvaluateSubtriangle(self,r0,r1,r2,z,n,Func):
        
        
        
        #lengths of all three sides of subdomain.
        l1 = np.subtract(r2,r1)
        l2 = np.subtract(r0,r2)
        l3 = np.subtract(r1,r0)
        cross1=fnp.fastCross(l1,l2)
        cross2=fnp.fastCross(l2,l3)
        cross3=fnp.fastCross(l3,l1)
        if fnp.fastNorm(cross1) == 0:
            if fnp.fastNorm(cross2) == 0:
                if fnp.fastNorm(cross3) == 0:
                    return 0
                else:
                    n_hat = cross3/(fnp.fastNorm(cross3))
            else:
                n_hat = cross2/(fnp.fastNorm(cross2))
        else:
            n_hat = cross1/(fnp.fastNorm(cross1)) #normal of subdomain
        #print(np.linalg.norm(np.cross(l1,l2)))
        
        Ap = np.dot(n_hat,cross1)/2 #area of subdomain
        h1 = np.multiply(2*Ap/(fnp.fastNorm(l1)**2),fnp.fastCross(l1,n_hat)) #height of subdomain. This will be along the y in our new coordinate system.  
        
        
        
        Gauss = np.polynomial.legendre.leggauss(self.order_duffy) #obtain weights and values for our gauss legendre polynomial of the given order
        #As we want gauss legendre from the intrval [0,1] instead of [-1,1] we need to transform the coordinates
        weight_new = np.multiply(0.5,Gauss[1]) 
        zeta_new = np.multiply(0.5,np.add(Gauss[0],1))
        length = len(zeta_new) #length of gaussian quadrature
        
        
        h1_hat = h1/fnp.fastNorm(h1) #direction of the heigth, this is the y_prime direction in our new coordinate system
        
        
        #upper and lower x values for the domain. As this has to of course be in our T2 triangle
        x_low = np.multiply(np.dot(n_hat,fnp.fastCross(h1_hat,l2)),np.subtract(1,zeta_new))
        x_up = -np.multiply(np.dot(n_hat,fnp.fastCross(h1_hat,l3)),np.subtract(1,zeta_new))        
        yp = np.multiply(fnp.fastNorm(h1),np.subtract(1,zeta_new)) #the y_prime calculation, this is the value thus in the h_hat direction
        
        x_hat = fnp.fastCross(h1_hat,n_hat) #determine the x_direction in our new coordinate system. After this point we have complete x_prime, y_prime, z_prime. Where the triangle lies in the x_prime, y_prime plane.
        
        #lower and upperbound of our substitution to cancle out the singularity.
        U_low = np.arcsinh(np.divide(x_low,np.sqrt(np.add(yp**2,z**2))))
        U_up = np.arcsinh(np.divide(x_up,np.sqrt(np.add(yp**2,z**2))))
        #print(r0)
        
        Int = 0
        num_j = np.complex128(0+1j) #again define the imaginairy unit
        
        #double summation over the triangle, for every y_prime there will be length points in the triangle. 
        for i in range(length):
            for j in range(length):
                U = U_low[j]*(1-zeta_new[i])+U_up[j]*zeta_new[i] #U points for triangle domain
                R = np.sqrt(yp[j]**2+z**2)*np.cosh(U) #This is the R equation to cancle singularity. See method duffy for reference.
                
                #x_test[i][j] = np.sqrt(R[i][j]**2-yp[j]**2)
                x_p = np.sqrt(yp[j]**2+z**2)*np.sinh(U) #x_prime location. needed to determine r_prime.
                rp = np.add(np.add(r0,np.multiply(x_p,x_hat)),np.multiply(yp[j],h1_hat)) #determine r_prime to be able to know where we are in actual coordinates instead of our own coordinate system
                const =  weight_new[i]*weight_new[j]*fnp.fastNorm(h1)*(U_up[j]-U_low[j]) #constant term see method duffy for reference
                Int += np.dot(n,n_hat)*const*Func(rp)*1/(4*np.pi)*np.exp(-num_j*self.k*R) #multiply everything and with the scalar green function without its singularity. 
        return Int
    
    
    #Function that calculates the integral over T2 with Duffy and Integral T1 with Dunavant. Every Dunavant point is an observation point for the dunavant method on T2.
    def int_bastest(self,T1,T2):
        
        Integrate = Basic_Integrator.Basic_Integrator(self.order_dunavant) #define integrator
        
        l1 = np.subtract(T2[1],T2[0])   #length of side 1 from vertex 2 to 3
        l2 = np.subtract(T2[2],T2[1])   #length of side 2 from vertex 1 to 2
        cross1=fnp.fastCross(l1,l2)
        n = np.divide(cross1,fnp.fastNorm(cross1)) #normal of triangle 2
        num_j = np.complex128(0+1j) #imaginairy unit
        
        length1 = fnp.fastNorm(np.subtract(T1[1],T1[0])) #length of common edge RWG
        length2 = fnp.fastNorm(np.subtract(T2[1],T2[0])) #length of common edge RWG
        A2 = Integrate.area_triangle(T2) #area of triangle 2
        A = Integrate.area_triangle(T1) #area of triangle 1
        
        
        basis = lambda r,rp: length2/(2*A2)*(np.subtract(rp,T2[2])) #define basis function
        Test = lambda r,rp: length1/(2*A)*(np.subtract(r,T1[2])) #define test function
        
        Div_basis_test = lambda r,rp: 1/(num_j*self.k)*length2/A2*length1/A #Divergence of test function, scalar so no need to make it a function
        basis_test = lambda r,rp: num_j*self.k*np.dot(basis(r,rp),Test(r,rp)) #function for the dot product of the tes and basis function. Easy way of implementing when providing both in first triangle, r will not matter as second triangle is only dependend on rp.
        

        proj = lambda r: np.add(r,np.multiply(np.divide(np.dot(n,np.subtract(T2[0],r)),fnp.fastNorm(n)**2),n)) #Function to calculate the location of the obervation point projected on triangle 2.
        z = lambda r: fnp.fastNorm(np.subtract(r,proj(r))) #function to calculate the height of observation point above projection point on T2.
        
        
        #proj = lambda r: np.subtract(np.subtract(r,T2[0]), np.dot(np.divide(np.dot(n,np.subtract(r,T2[0])),np.dot(n,n)),n))
        
        #Function to integrate the first integral this is the hypersingular equation where the divergences are taken into account. Therefroe the function given is 1 as the divergence is multiplied later as this is a constant.
        #The summation is implemented three times to calculate all the three domains of triangle 2 using duffy.
        Int_triangle = lambda r: self.EvaluateSubtriangle(proj(r),T2[1],T2[2],z(r),n,lambda rp: Div_basis_test(r,rp)) + self.EvaluateSubtriangle(proj(r),T2[2],T2[0],z(r),n,lambda rp: Div_basis_test(r,rp)) + self.EvaluateSubtriangle(proj(r),T2[0],T2[1],z(r),n,lambda rp: Div_basis_test(r,rp))
        I_HS = Integrate.int_triangle(Int_triangle,T1)
        
        
        #Function to integrate the second integral the Singular integral, here the basis and test function are given. This needs to be done three times for all the three triangular domains for duffy.
        Int_triangle2 = lambda r: self.EvaluateSubtriangle(proj(r),T2[1],T2[2],z(r),n,lambda rp: basis_test(r,rp)) + self.EvaluateSubtriangle(proj(r),T2[2],T2[0],z(r),n,lambda rp: basis_test(r,rp)) + self.EvaluateSubtriangle(proj(r),T2[0],T2[1],z(r),n,lambda rp: basis_test(r,rp))
        I_S = Integrate.int_triangle(Int_triangle2,T1)
            
        
        return I_HS+I_S
