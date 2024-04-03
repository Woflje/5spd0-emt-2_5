import numpy as np
class Integrator_base_test_duffy():
    def __init__(self,k,Gauss_order,dunavant_weight):
        self.k = k
        self.dunavant_weight = dunavant_weight
        self.dunavant_weight_expanded = self.dunavant_weight[:,np.newaxis]
        self.dunavant_length = len(self.dunavant_weight)

        Gauss = np.polynomial.legendre.leggauss(Gauss_order)
        self.Gauss_weight = 0.5*Gauss[1]
        self.Gauss_zeta = 0.5*(Gauss[0]+1)
        self.Gauss_zeta_complement = 1-self.Gauss_zeta
        self.Gauss_weight_outer = np.outer(self.Gauss_weight, self.Gauss_weight)

    def EvaluateSubtriangle(self, r0, r1, r2, z, n, integrand):
        l1 = r2 - r1
        l2 = r0 - r2
        l3 = r1 - r0

        cross_l1_l2 = np.cross(l1[None, :], l2)
        cross_l2_l3 = np.cross(l2, l3)
        cross_l3_l1 = np.cross(l3, l1[None, :])

        cross_products = np.stack([cross_l1_l2, cross_l2_l3, cross_l3_l1], axis=0)
        norms = np.linalg.norm(cross_products, axis=2)

        # Identify colinear points by checking if all norms per row of r0 are zero
        colinear_indices = np.where(np.all(norms == 0, axis=0))[0]

        n_hat = np.zeros((self.dunavant_length, 3))

        # Iterate over each set of cross products and their corresponding norms
        for i in range(3):
            valid = norms[i] != 0
            if np.any(valid):
                # Normalize only where norms are not zero
                n_hat[valid] = cross_products[i, valid] / norms[i, valid][:, np.newaxis]
        # Pre-vectorization, this function would return 0 if no non-zero norms were found for a position
        # As this is not possible with a vector of positions, the results of the function are set to 0 at colinear_indices

        # Subdomain area
        A_sub = np.einsum('ij,ij->i', n_hat, cross_l1_l2) / 2

        # Subdomain height
        h1 = 2 * A_sub[:, None] / np.linalg.norm(l1)**2 * np.cross(l1, n_hat)
        h1_norm = np.linalg.norm(h1, axis=1)
        h1_hat = h1 / h1_norm[:, None]

        # y_prime
        yp = h1_norm[:, None] * self.Gauss_zeta_complement
        r_eff = np.sqrt(yp**2 + z**2)

        # Lower and upper bounds
        x_low = np.einsum('ij,ij->i', n_hat, np.cross(h1_hat, l2))[:, None] * self.Gauss_zeta_complement
        x_up = -np.einsum('ij,ij->i', n_hat, np.cross(h1_hat, l3))[:, None] * self.Gauss_zeta_complement

        U_low = np.arcsinh(x_low / r_eff)
        U_up = np.arcsinh(x_up / r_eff)

        # Duffy integration, as explained in the paper on Duffy (See references)
        U = np.einsum('i,jk->jik', 1 - self.Gauss_zeta, U_low) + np.einsum('i,jk->jik', self.Gauss_zeta, U_up)
        R = r_eff[:,np.newaxis,:] * np.cosh(U)
        n_dot_nhat = np.dot(n, n_hat.T)
        
        U_diff = U_up - U_low

        const = self.Gauss_weight_outer * h1_norm[:,np.newaxis,np.newaxis] * U_diff[:,np.newaxis,:]

        exp_term = np.exp(-1j * self.k * R)

        if not callable(integrand):
            # Multiply everything and with the scalar green function without its singularity. 
            integral_result = np.sum(n_dot_nhat[:, None, None] * const * integrand / (4 * np.pi) * exp_term, axis=(1, 2))
        else:
            x_hat = np.cross(h1_hat, n_hat)
            x_p = r_eff[:,np.newaxis,:]*np.sinh(U) # x_prime location. needed to determine r_prime.
            outer_product = (yp[:, :, np.newaxis] * h1_hat[:, np.newaxis, :]).reshape(yp.shape[0], 1, yp.shape[1], 3)
            m_x_p_x_hat = x_p[:, :, :, np.newaxis] * x_hat[:, np.newaxis, np.newaxis, :]
            rp = r0[:, None, None] + m_x_p_x_hat + outer_product

            # Multiply everything and with the scalar green function without its singularity. 
            integral_result = np.sum(n_dot_nhat[:, None, None] * const * integrand(rp) / (4 * np.pi) * exp_term, axis=(1, 2))

        integral_result[colinear_indices]=0
        return integral_result
    
    def integrate_base_test(self,T1,T2,A1,A2,dunavant_positions1):
        vector_2_3 = np.subtract(T2[1],T2[0])
        vector_1_2 = np.subtract(T2[2],T2[1])
        cross1=np.cross(vector_2_3,vector_1_2)
        n = np.divide(cross1,np.linalg.norm(cross1)) # Normal of triangle 2
        
        common_edge_length_t1 = np.linalg.norm(T1[1]-T1[0]) #length of common edge RWG
        common_edge_length_t2 = np.linalg.norm(vector_2_3) #length of common edge RWG
        
        Test = common_edge_length_t1/(2*A1)*(dunavant_positions1-T1[2])
        
        # Divergence of test function
        Divergence_basis_test = 1/(1j*self.k)*common_edge_length_t2/A2*common_edge_length_t1/A1 
        
        difference = T2[0]-dunavant_positions1

        # Calculate the dot product of n with each row of the difference.
        dot_product = np.dot(difference, n.reshape(-1, 1)).reshape(-1)

        # Location of the observation point projected on triangle 2
        proj = dunavant_positions1 + np.outer(dot_product / np.linalg.norm(n)**2, n)

        # Height of observation point above projection point on triangle 2.
        z = np.linalg.norm(dunavant_positions1-proj, axis=1)[:,None]

        def basis_test():
            def integrand(rp):
                basis = common_edge_length_t2 / (2*A2)*(rp-T2[2])
                return 1j*self.k*np.einsum('pqkc,pc->pqk',basis,Test)
            return integrand

        # Function to integrate the first integral this, is the hypersingular equation where the divergences are taken into account.
        # The summation is implemented three times to calculate all the three domains of triangle 2 using Duffy.

        # Hyper Singular
        Integrand_triangle_1 = sum([
            self.EvaluateSubtriangle(proj,T2[1],T2[2],z,n,Divergence_basis_test),
            self.EvaluateSubtriangle(proj,T2[2],T2[0],z,n,Divergence_basis_test),
            self.EvaluateSubtriangle(proj,T2[0],T2[1],z,n,Divergence_basis_test)
        ])

        integrand_function = basis_test()

        # Singular
        Integrand_triangle_2 = sum([
            self.EvaluateSubtriangle(proj,T2[1],T2[2],z,n,integrand_function),
            self.EvaluateSubtriangle(proj,T2[2],T2[0],z,n,integrand_function),
            self.EvaluateSubtriangle(proj,T2[0],T2[1],z,n,integrand_function)
        ])

        I_HS = np.sum(Integrand_triangle_1[:, None] * self.dunavant_weight_expanded)
        I_S = np.sum(Integrand_triangle_2[:, None] * self.dunavant_weight_expanded)
        result = A1*(I_HS+I_S)
        return result