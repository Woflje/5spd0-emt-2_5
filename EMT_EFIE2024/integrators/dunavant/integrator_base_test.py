import numpy as np

class Integrator_base_test_dunavant():
    def __init__(self,k,dunavant_weight):
        self.k = k
        self.dunavant_weight = dunavant_weight
    
    def _scalar_green(self, r, rp):
        # Calculate the norm between each pair of r and rp points
        norm = np.linalg.norm(r[:, None, :] - rp[None, :, :], axis=2)
        return np.exp(-1j * self.k * norm) / (4 * np.pi * norm)
    
    def integrate_base_test(self,T1,T2,A1,A2,dunavant_pos1,dunavant_pos2):
        length1 = np.linalg.norm(T1[1]-T1[0]) # Length common edge rwg 1
        length2 = np.linalg.norm(T2[1]-T2[0]) # Length common edge rwg 2
        
        # Precompute scalar green's function for all pairs of points
        green_matrix1 = self._scalar_green(dunavant_pos1, dunavant_pos2)

        # Dunavant integral method
        integrals_part1 = 1/(1j*self.k)*(length1*length2/(A1*A2))*green_matrix1
        result_integral1 = A1*np.sum(np.sum(integrals_part1*self.dunavant_weight*A2,axis=1)*self.dunavant_weight)

        term1 = length1/(2*A1)*(dunavant_pos1-T1[2])
        term2 = length2/(2*A2)*(dunavant_pos2-T2[2])

        dot_products = np.einsum('ij,kj->ik',term1,term2)
        integrals_part2 = 1j*self.k*dot_products*green_matrix1
        result_integral2 = A1*np.sum(np.sum(integrals_part2*self.dunavant_weight*A2,axis=1)*self.dunavant_weight)

        result = result_integral1+result_integral2

        return result