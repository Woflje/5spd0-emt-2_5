import numpy as np
import FastNP as fnp
from loose import calculate_triangle_area
from Integration import int_triangle, int_triangle_duffy, int_triangle_duffy_scalar


class integrator_bastest():
    def __init__(self,dunavant_integrator,k,order_duffy,order_dunavant,dunavant_values,dunavant_weight):
        print("Initialize integration Duffy method")
        self.dunavant_integrator = dunavant_integrator
        self.k = k # wave vector
        self.order_duffy = order_duffy #oder of the duffy method
        self.order_dunavant = order_dunavant #order of the dunvant method
        self.dunavant_values = dunavant_values
        self.dunavant_weight = dunavant_weight
        self.dunavant_weight_expanded = self.dunavant_weight[:,np.newaxis]

        self.Gauss = np.polynomial.legendre.leggauss(self.order_duffy)
        self.weight_new = 0.5*self.Gauss[1]
        self.zeta_new = 0.5*(self.Gauss[0]+1)
        self.zeta_new_ = 1-self.zeta_new
        self.weight_outer = np.outer(self.weight_new, self.weight_new)

    def EvaluateSubtriangle_loop(self,r0,r1,r2,z,n,func):
        #r0 (P,3) array of 3rd vertex positions
        #r1 (3,) constant vertex 1
        #r2 (3,) constant vertex 2
        #r1 and r2 make up the constant side of the triangles to evaluate
        #z (P,)
        #n (3,) normal
        #func takes a single value but I vectorized it so it can take (P,)
        #return value of this function is (P,)
        P = len(r0)
        Int=[]
        for r in range(P):
        # lengths of sides
            l1,l2,l3=r2-r1, r0[r]-r2, r1-r0[r]
            # cross products
            cross = [np.cross(l1,l2), np.cross(l2,l3), np.cross(l3,l1)]
            norms = [np.linalg.norm(x) for x in cross]
            
            # find the first non-zero normal vector
            n_hat = next((x/norm for x, norm in zip(cross, norms) if norm != 0), None)
            if n_hat is None:
                Int.append(0)
                continue
            
            # Subdomain area
            A_sub = np.dot(n_hat, cross[0])/2
            
            # Subdomain height
            h1 = 2*A_sub/np.square(np.linalg.norm(l1))*np.cross(l1, n_hat)
            h1_norm = np.linalg.norm(h1)
            h1_hat = h1/h1_norm

            Q = 5

            # Prepare yp (y_prime)
            yp = h1_norm*self.zeta_new_

            # Prepare effective radius
            r_eff = np.sqrt(yp**2+z[r]**2) #shape (5,)

            # Lower and upper bounds
            x_low = np.dot(n_hat, np.cross(h1_hat,l2))*self.zeta_new_
            x_up = -np.dot(n_hat, np.cross(h1_hat,l3))*self.zeta_new_
            U_low = np.arcsinh(x_low/r_eff)
            U_up = np.arcsinh(x_up/r_eff)

            # Double summation
            U = np.outer(self.zeta_new_, U_low) + np.outer(self.zeta_new, U_up) #shape (5,5)
            R = r_eff*np.cosh(U) 
            
            # Integration constant
            const = self.weight_outer*h1_norm*(U_up-U_low)

            # Final integration
            if not callable(func):
                Int_temp = np.sum((np.dot(n,n_hat)*const*func/(4*np.pi)*np.exp(-1j*self.k*R)))
            else:
                # x_hat direction
                x_hat = np.cross(h1_hat,n_hat)
                x_p = r_eff*np.sinh(U)
                # rp and dimension struggle
                r0_expanded = np.expand_dims(r0[r], axis=0)
                rp = r0_expanded + np.tensordot(x_p, x_hat, axes=0) + np.outer(yp, h1_hat).reshape(-1,1,3)
                Int_temp = np.sum((np.dot(n,n_hat)*const*func(rp)/(4*np.pi)*np.exp(-1j*self.k*R)))
            Int.append(Int_temp)
        return np.array(Int)

    def EvaluateSubtriangle(self,r0,r1,r2,z,n,func):
        # lengths of sides
        l1,l2,l3=r2-r1, r0-r2, r1-r0
        # cross products
        cross = [np.cross(l1,l2), np.cross(l2,l3), np.cross(l3,l1)]
        norms = [np.linalg.norm(x) for x in cross]
        
        # find the first non-zero normal vector
        n_hat = next((x/norm for x, norm in zip(cross, norms) if norm != 0), None)
        if n_hat is None:
            return 0 # Sides are collinear, no valid subtriangle
        
        # Subdomain area
        A_sub = np.dot(n_hat, cross[0])/2
        
        # Subdomain height
        h1 = 2*A_sub/np.square(np.linalg.norm(l1))*np.cross(l1, n_hat)
        h1_norm = np.linalg.norm(h1)
        h1_hat = h1/h1_norm

        # Prepare yp (y_prime)
        yp = h1_norm*self.zeta_new_

        # Prepare effective radius
        r_eff = np.sqrt(yp**2+z**2)

        # Lower and upper bounds
        x_low = np.dot(n_hat, np.cross(h1_hat,l2))*self.zeta_new_
        x_up = -np.dot(n_hat, np.cross(h1_hat,l3))*self.zeta_new_
        U_low = np.arcsinh(x_low/r_eff)
        U_up = np.arcsinh(x_up/r_eff)

        # Double summation
        U = np.outer(self.zeta_new_, U_low) + np.outer(self.zeta_new, U_up)
        R = r_eff*np.cosh(U)
        
        # Integration constant
        const = self.weight_outer*h1_norm*(U_up-U_low)

        # Final integration
        if not callable(func):
            Int = np.sum((np.dot(n,n_hat)*const*func/(4*np.pi)*np.exp(-1j*self.k*R)))
        else:

            # x_hat direction
            x_hat = np.cross(h1_hat,n_hat)
            x_p = r_eff*np.sinh(U)
            # rp and dimension struggle
            r0_expanded = np.expand_dims(r0, axis=0)
            rp = r0_expanded + np.tensordot(x_p, x_hat, axes=0) + np.outer(yp, h1_hat).reshape(-1,1,3)
            Int = np.sum((np.dot(n,n_hat)*const*func(rp)/(4*np.pi)*np.exp(-1j*self.k*R)))
        return Int

        #Function that computes the integral over one subdomain. If called three times the whole triangle is integrated using duffy
    def EvaluateSubtriangle_even_newer(self, r0, r1, r2, z, n, integrand):
        P = r0.shape[0]
        # Basic vector operations
        l1 = r2 - r1
        l2 = r0 - r2
        l3 = r1 - r0

        # Vectorized cross products
        cross_l1_l2 = np.cross(l1[None, :], l2)
        cross_l2_l3 = np.cross(l2, l3)
        cross_l3_l1 = np.cross(l3, l1[None, :])

        cross_products = np.stack([cross_l1_l2, cross_l2_l3, cross_l3_l1], axis=0)  # Shape: (3, P, 3)
        norms = np.linalg.norm(cross_products, axis=2)  # Shape: (3, P)

        # Identify colinear points by checking if all norms are zero
        colinear_indices = np.where(np.all(norms == 0, axis=0))[0]

        # Initialize n_hat with zeros with the correct shape (P, 3)
        n_hat = np.zeros((P, 3))

        # Iterate over each set of cross products and their corresponding norms
        for i in range(3):
            valid = norms[i] != 0
            if np.any(valid):
                # Normalize only where norms are not zero
                n_hat[valid] = cross_products[i, valid] / norms[i, valid][:, np.newaxis]

        # Subdomain area
        A_sub = np.einsum('ij,ij->i', n_hat, cross_l1_l2) / 2

        # Subdomain height
        h1 = 2 * A_sub[:, None] / np.linalg.norm(l1)**2 * np.cross(l1, n_hat)
        h1_norm = np.linalg.norm(h1, axis=1)
        h1_hat = h1 / h1_norm[:, None]

        # Gauss-Legendre Quadrature for integration points and weights
        Q = 5

        # Prepare yp (y_prime) and effective radius in a vectorized way
        yp = h1_norm[:, None] * self.zeta_new_
        r_eff = np.sqrt(yp**2 + z**2)

        # Lower and upper bounds
        x_low = np.einsum('ij,ij->i', n_hat, np.cross(h1_hat, l2))[:, None] * self.zeta_new_
        x_up = -np.einsum('ij,ij->i', n_hat, np.cross(h1_hat, l3))[:, None] * self.zeta_new_

        U_low = np.arcsinh(x_low / r_eff)
        U_up = np.arcsinh(x_up / r_eff)
        
        # Double summation and Final integration
        U = np.empty((P, Q, Q))
        for i in range(P):
            U[i] = np.outer(self.zeta_new_, U_low[i]) + np.outer(self.zeta_new, U_up[i])

        # Correcting the calculation for R with the right shape
        R = r_eff[:, None] * np.cosh(U)

        # Integration constant
        n_dot_n_hat = np.sum(n * n_hat, axis=1)  # Shape: (P,)
        
        # Integration constant with corrected shapes
        h1_norm_expanded = h1_norm[:, None, None]  # Expand h1_norm to shape (P, 1, 1) for broadcasting
        U_diff = U_up - U_low  # This subtraction is fine, both have shape (P, Q)

        # Now, we need to broadcast U_diff correctly to a (P, Q, Q) shape. The tricky part is doing this in a way that matches with weight_outer
        # We expand U_diff for broadcasting across the Q dimension not covered by the subtraction itself
        U_diff_expanded = U_diff[:, :, None]  # Now has shape (P, Q, 1), ready to broadcast across the second Q of weight_outer
        # Finally, multiply everything together, ensuring correct broadcasting
        const = self.weight_outer * h1_norm_expanded * U_diff_expanded
        
        exp_term = np.exp(-1j * self.k * R)

        # Integration over the Gauss points
        if not callable(func):
            Int = np.sum(n_dot_n_hat[:, None, None] * const * func / (4 * np.pi) * exp_term, axis=(1, 2))
        else:
            # Variable input for callable func
            x_hat = np.cross(h1_hat, n_hat)  # Shape: (P, 3)
            x_p = r_eff[:, :, None] * np.sinh(U)  # Shape: (P, Q, Q), needs reshaping for correct application
            
            rp = []
            for i in range(P):
                r0_expanded = np.expand_dims(r0[i], axis=0)
                rp.append(r0_expanded + np.tensordot(x_p[i], x_hat[i], axes=0) + np.outer(yp[i], h1_hat[i]).reshape(-1,1,3))
            rp = np.array(rp)
            Int = np.sum((np.dot(n, n_hat.T)[:, None, None] * const * func(rp) / (4 * np.pi) * np.exp(-1j * self.k * R)), axis=(1, 2))
        Int[colinear_indices]=0
        return Int
    def EvaluateSubtriangle_new(self, r0, r1, r2, z, n, func):
        P = r0.shape[0]
        # Basic vector operations
        l1 = r2 - r1
        l2 = r0 - r2
        l3 = r1 - r0

        # Vectorized cross products
        cross_l1_l2 = np.cross(l1[None, :], l2)
        cross_l2_l3 = np.cross(l2, l3)
        cross_l3_l1 = np.cross(l3, l1[None, :])

        cross_products = np.stack([cross_l1_l2, cross_l2_l3, cross_l3_l1], axis=0)  # Shape: (3, P, 3)
        norms = np.linalg.norm(cross_products, axis=2)  # Shape: (3, P)

        # Identify colinear points by checking if all norms are zero
        colinear_indices = np.where(np.all(norms == 0, axis=0))[0]

        # Initialize n_hat with zeros with the correct shape (P, 3)
        n_hat = np.zeros((P, 3))

        # Iterate over each set of cross products and their corresponding norms
        for i in range(3):
            valid = norms[i] != 0
            if np.any(valid):
                # Normalize only where norms are not zero
                n_hat[valid] = cross_products[i, valid] / norms[i, valid][:, np.newaxis]

        # Subdomain area
        A_sub = np.einsum('ij,ij->i', n_hat, cross_l1_l2) / 2

        # Subdomain height
        h1 = 2 * A_sub[:, None] / np.linalg.norm(l1)**2 * np.cross(l1, n_hat)
        h1_norm = np.linalg.norm(h1, axis=1)
        h1_hat = h1 / h1_norm[:, None]

        # Gauss-Legendre Quadrature for integration points and weights
        Q = 5

        # Prepare yp (y_prime) and effective radius in a vectorized way
        yp = h1_norm[:, None] * self.zeta_new_
        r_eff = np.sqrt(yp**2 + z**2)

        # Lower and upper bounds
        x_low = np.einsum('ij,ij->i', n_hat, np.cross(h1_hat, l2))[:, None] * self.zeta_new_
        x_up = -np.einsum('ij,ij->i', n_hat, np.cross(h1_hat, l3))[:, None] * self.zeta_new_

        U_low = np.arcsinh(x_low / r_eff)
        U_up = np.arcsinh(x_up / r_eff)
        
        # Double summation and Final integration
        # U = np.empty((P, Q, Q))
        # for i in range(P):
        #     U[i] = np.outer(1-self.zeta_new, U_low[i]) + np.outer(self.zeta_new, U_up[i])

        U = np.einsum('i,jk->jik', 1 - self.zeta_new, U_low) + np.einsum('i,jk->jik', self.zeta_new, U_up)

        # Correcting the calculation for R with the right shape
        R = r_eff[:, None] * np.cosh(U)

        # Integration constant
        n_dot_n_hat = np.sum(n * n_hat, axis=1)  # Shape: (P,)
        
        # Integration constant with corrected shapes
        h1_norm_expanded = h1_norm[:, None, None]  # Expand h1_norm to shape (P, 1, 1) for broadcasting
        U_diff = U_up - U_low  # This subtraction is fine, both have shape (P, Q)

        # Now, we need to broadcast U_diff correctly to a (P, Q, Q) shape. The tricky part is doing this in a way that matches with weight_outer
        # We expand U_diff for broadcasting across the Q dimension not covered by the subtraction itself
        U_diff_expanded = U_diff[:, :, None]  # Now has shape (P, Q, 1), ready to broadcast across the second Q of weight_outer
        # Finally, multiply everything together, ensuring correct broadcasting
        const = self.weight_outer * h1_norm_expanded * U_diff_expanded
        
        exp_term = np.exp(-1j * self.k * R)

        # Integration over the Gauss points
        if not callable(func):
            Int = np.sum(n_dot_n_hat[:, None, None] * const * func / (4 * np.pi) * exp_term, axis=(1, 2))
        else:
            x_hat = np.cross(h1_hat, n_hat)
            x_p = r_eff[..., None] * np.sinh(U)
            outer_product = yp[..., None] * h1_hat[:, None]
            rp = r0[:, None, None] + x_p[..., None] * x_hat[:, None, None] + outer_product[:, :, None]
            n_dot_nhat = np.dot(n, n_hat.T)
            Int = np.sum(n_dot_nhat[:, None, None] * const * func(rp) / (4 * np.pi) * exp_term, axis=(1, 2))
        Int[colinear_indices]=0
        return Int
    
    def int_bastest_new(self,T1,T2,A1,A2,dunavant_positions1):
        vector_2_3 = np.subtract(T2[1],T2[0])   #length of side 1 from vertex 2 to 3
        vector_1_2 = np.subtract(T2[2],T2[1])   #length of side 2 from vertex 1 to 2
        cross1=np.cross(vector_2_3,vector_1_2)
        n = np.divide(cross1,np.linalg.norm(cross1)) #normal of triangle 2
        
        common_edge_length_t1 = np.linalg.norm(np.subtract(T1[1],T1[0])) #length of common edge RWG
        common_edge_legnth_t2 = np.linalg.norm(vector_2_3) #length of common edge RWG\
        
        Test = common_edge_length_t1/(2*A1)*(np.subtract(dunavant_positions1,T1[2]))

        basis = lambda rp: common_edge_legnth_t2/(2*A2)*(np.subtract(rp,T2[2])) #define basis function
        
        Div_basis_test = 1/(1j*self.k)*common_edge_legnth_t2/A2*common_edge_length_t1/A1 #Divergence of test function, scalar so no need to make it a function
        basis_test = lambda rp: 1j*self.k*np.einsum('pqkc,pc->pqk', basis(rp), Test) #function for the dot product of the tes and basis function. Easy way of implementing when providing both in first triangle, r will not matter as second triangle is only dependend on rp.
        
        # location of the obervation point projected on triangle 2.
        difference = np.subtract(T2[0], dunavant_positions1) # This will be of shape (P, 3)

        # Calculate the dot product of n with each row of the difference. The reshape ensures compatibility.
        dot_product = np.dot(difference, n.reshape(-1, 1)).reshape(-1) # Reshape to (P,) for broadcasting

        # Divide by the squared norm of n, which is a scalar
        division = dot_product / np.linalg.norm(n)**2

        # Multiply by n (broadcasted automatically) and add to the original positions
        proj = dunavant_positions1 + np.outer(division, n)
        z = np.linalg.norm(np.subtract(dunavant_positions1, proj), axis=1)[:,None]
        
        #Function to integrate the first integral this is the hypersingular equation where the divergences are taken into account. Therefroe the function given is 1 as the divergence is multiplied later as this is a constant.
        #The summation is implemented three times to calculate all the three domains of triangle 2 using duffy.

        Int_triangle1 = self.EvaluateSubtriangle_new(proj,T2[1],T2[2],z,n,Div_basis_test)\
            + self.EvaluateSubtriangle_new(proj,T2[2],T2[0],z,n,Div_basis_test)\
            + self.EvaluateSubtriangle_new(proj,T2[0],T2[1],z,n,Div_basis_test)
        # I_HS = int_triangle_duffy_scalar(Int_triangle1,A1,self.dunavant_weight)
        I_HS = sum(np.multiply(Int_triangle1[:,np.newaxis],self.dunavant_weight_expanded))[0]
        #Function to integrate the second integral the Singular integral, here the basis and test function are given. This needs to be done three times for all the three triangular domains for duffy.

        Int_triangle2 = lambda rp: self.EvaluateSubtriangle_new(proj,T2[1],T2[2],z,n,basis_test)\
            + self.EvaluateSubtriangle_new(proj,T2[2],T2[0],z,n,basis_test)\
            + self.EvaluateSubtriangle_new(proj,T2[0],T2[1],z,n,basis_test)
        # I_S = int_triangle_duffy(Int_triangle2,A1,dunavant_positions1,self.dunavant_weight)
        I_S = sum(np.multiply(Int_triangle2(dunavant_positions1)[:,np.newaxis],self.dunavant_weight_expanded))[0]
            
        return A1*(I_HS+I_S)
    
    #Function that calculates the integral over T2 with Duffy and Integral T1 with Dunavant. Every Dunavant point is an observation point for the dunavant method on T2.
    def int_bastest(self,T1,T2):
        l1 = np.subtract(T2[1],T2[0])   #length of side 1 from vertex 2 to 3
        l2 = np.subtract(T2[2],T2[1])   #length of side 2 from vertex 1 to 2
        cross1=np.cross(l1,l2)
        n = np.divide(cross1,np.linalg.norm(cross1)) #normal of triangle 2
        num_j = np.complex128(0+1j) #imaginairy unit
        
        length1 = np.linalg.norm(np.subtract(T1[1],T1[0])) #length of common edge RWG
        length2 = np.linalg.norm(np.subtract(T2[1],T2[0])) #length of common edge RWG\
        A1 = calculate_triangle_area(T1) #area of triangle 1
        A2 = calculate_triangle_area(T2) #area of triangle 2
        
        basis = lambda r,rp: length2/(2*A2)*(np.subtract(rp,T2[2])) #define basis function
        Test = lambda r,rp: length1/(2*A1)*(np.subtract(r,T1[2])) #define test function
        
        Div_basis_test = 1/(num_j*self.k)*length2/A2*length1/A1 #Divergence of test function, scalar so no need to make it a function
        basis_test = lambda r,rp: num_j*self.k*np.dot(basis(r,rp),Test(r,rp)) #function for the dot product of the tes and basis function. Easy way of implementing when providing both in first triangle, r will not matter as second triangle is only dependend on rp.
        

        proj = lambda r: np.add(r,np.multiply(np.divide(np.dot(n,np.subtract(T2[0],r)),np.linalg.norm(n)**2),n)) #Function to calculate the location of the obervation point projected on triangle 2.
        z = lambda r: np.linalg.norm(np.subtract(r,proj(r))) #function to calculate the height of observation point above projection point on T2.
        
        #Function to integrate the first integral this is the hypersingular equation where the divergences are taken into account. Therefroe the function given is 1 as the divergence is multiplied later as this is a constant.
        #The summation is implemented three times to calculate all the three domains of triangle 2 using duffy.

        dunavant_positions_T1=np.dot(self.dunavant_values,T1)

        Int_triangle1 = lambda r: self.EvaluateSubtriangle(proj(r),T2[1],T2[2],z(r),n,Div_basis_test) + self.EvaluateSubtriangle(proj(r),T2[2],T2[0],z(r),n,Div_basis_test) + self.EvaluateSubtriangle(proj(r),T2[0],T2[1],z(r),n,Div_basis_test)
        I_HS = int_triangle(Int_triangle1,A1,dunavant_positions_T1,self.dunavant_weight)
        
        
        #Function to integrate the second integral the Singular integral, here the basis and test function are given. This needs to be done three times for all the three triangular domains for duffy.
        Int_triangle2 = lambda r: self.EvaluateSubtriangle(proj(r),T2[1],T2[2],z(r),n,lambda rp: basis_test(r,rp)) + self.EvaluateSubtriangle(proj(r),T2[2],T2[0],z(r),n,lambda rp: basis_test(r,rp)) + self.EvaluateSubtriangle(proj(r),T2[0],T2[1],z(r),n,lambda rp: basis_test(r,rp))
        I_S = int_triangle(Int_triangle2,A1,dunavant_positions_T1,self.dunavant_weight)
            
        
        return I_HS+I_S
