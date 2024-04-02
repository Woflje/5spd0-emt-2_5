import numpy as np
from Timer import Timer
class Integrator_base_test_dunavant():
    def __init__(self,k,dunavant_weight):
        self.k = k # wave vector
        self.dunavant_weight = dunavant_weight
    
    def _scalar_green(self, r, rp):
        # Calculate the norm between each pair of r and rp points
        norm = np.linalg.norm(r[:, None, :] - rp[None, :, :], axis=2)
        return np.exp(-1j * self.k * norm) / (4 * np.pi * norm)
    
    def integrate_base_test(self,T1,T2,A1,A2,dunavant_pos1,dunavant_pos2): #initialize the basic integrators
        length1 = np.linalg.norm(T1[1]-T1[0]) #length common edge rwg 1
        length2 = np.linalg.norm(T2[1]-T2[0]) #length common edge rwg 2
        
        # Precompute scalar green's function for all pairs of points
        green_matrix1 = self._scalar_green(dunavant_pos1, dunavant_pos2)

        integrals_part1 = 1/(1j*self.k)*(length1*length2/(A1*A2))*green_matrix1
        result_integral1 = A1*np.sum(np.sum(integrals_part1*self.dunavant_weight*A2,axis=1)*self.dunavant_weight)

        term1 = length1/(2*A1)*(dunavant_pos1-T1[2])
        term2 = length2/(2*A2)*(dunavant_pos2-T2[2])

        dot_products = np.einsum('ij,kj->ik',term1,term2)
        integrals_part2 = 1j*self.k*dot_products*green_matrix1
        result_integral2 = A1*np.sum(np.sum(integrals_part2*self.dunavant_weight*A2,axis=1)*self.dunavant_weight)

        result = result_integral1+result_integral2

        return result
    
    def integrate_base_test_loop(self,T1,T2,A1,A2,dunavant_pos1,dunavant_pos2):
        length1 = np.linalg.norm(T1[1]-T1[0]) #length common edge rwg 1
        length2 = np.linalg.norm(T2[1]-T2[0]) #length common edge rwg 2

        temp1 = []
        #part 1
        for r in dunavant_pos1:
            temp2 = []
            for rp in dunavant_pos2:
                norm = np.linalg.norm(np.subtract(r,rp))
                scalar_green = np.exp(-1j*self.k*norm)/(4*np.pi*norm)
                part1 = 1/(1j*self.k)*length2/A2*length1/A1*scalar_green
                temp2.append(part1)
            temp1.append(A2*sum(np.multiply(temp2,self.dunavant_weight)))
        result_integral1 = A1*sum(np.multiply(temp1,self.dunavant_weight))

        temp1 = []
        #part 2
        for r in dunavant_pos1:
            temp2 = []
            for rp in dunavant_pos2:
                norm = np.linalg.norm(np.subtract(r,rp))
                scalar_green = np.exp(-1j*self.k*norm)/(4*np.pi*norm)
                term1 = length1/(2*A1)*(r-T1[2])
                term2 = length1/(2*A2)*(rp-T2[2])
                part2 = 1j*self.k*np.dot(term1,term2)*scalar_green
                temp2.append(part2)
            temp1.append(A2*sum(np.multiply(temp2,self.dunavant_weight)))
        result_integral2 = A1*sum(np.multiply(temp1,self.dunavant_weight))
        
        result = result_integral1+result_integral2
        return result