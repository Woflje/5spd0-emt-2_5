import numpy as np

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
        self.dunavant_length = len(self.dunavant_weight)

        self.Gauss = np.polynomial.legendre.leggauss(self.order_duffy)
        self.Gauss_weight = 0.5*self.Gauss[1]
        self.Gauss_zeta = 0.5*(self.Gauss[0]+1)
        self.Gauss_zeta_complement = 1-self.Gauss_zeta
        self.Gauss_weight_outer = np.outer(self.Gauss_weight, self.Gauss_weight)
        self.Gauss_length = len(self.Gauss_zeta)

    def EvaluateSubtriangle_new(self, r0, r1, r2, z, n, integrand):
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
        n_hat = np.zeros((self.dunavant_length, 3))

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

        # Prepare yp (y_prime) and effective radius in a vectorized way
        yp = h1_norm[:, None] * self.Gauss_zeta_complement
        r_eff = np.sqrt(yp**2 + z**2)

        # Lower and upper bounds
        x_low = np.einsum('ij,ij->i', n_hat, np.cross(h1_hat, l2))[:, None] * self.Gauss_zeta_complement
        x_up = -np.einsum('ij,ij->i', n_hat, np.cross(h1_hat, l3))[:, None] * self.Gauss_zeta_complement

        U_low = np.arcsinh(x_low / r_eff)
        U_up = np.arcsinh(x_up / r_eff)

        U = np.einsum('i,jk->jik', 1 - self.Gauss_zeta, U_low) + np.einsum('i,jk->jik', self.Gauss_zeta, U_up)

        R = r_eff[:, None] * np.cosh(U)

        n_dot_n_hat = np.sum(n * n_hat, axis=1)
        
        h1_norm_expanded = h1_norm[:, None, None]
        U_diff = U_up - U_low

        const = self.Gauss_weight_outer * h1_norm_expanded * U_diff[:, :, None]
        
        exp_term = np.exp(-1j * self.k * R)

        # Integration over the Gauss points
        if not callable(integrand):
            Int = np.sum(n_dot_n_hat[:, None, None] * const * integrand / (4 * np.pi) * exp_term, axis=(1, 2))
        else:
            x_hat = np.cross(h1_hat, n_hat)
            x_p = r_eff[..., None] * np.sinh(U)
            outer_product = yp[..., None] * h1_hat[:, None]
            rp = r0[:, None, None] + x_p[..., None] * x_hat[:, None, None] + outer_product[:, :, None]
            n_dot_nhat = np.dot(n, n_hat.T)
            Int = np.sum(n_dot_nhat[:, None, None] * const * integrand(rp) / (4 * np.pi) * exp_term, axis=(1, 2))
        Int[colinear_indices]=0
        return Int
    
    def int_bastest(self,T1,T2,A1,A2,dunavant_positions1):
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