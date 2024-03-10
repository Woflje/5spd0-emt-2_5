import numpy as np
from Dunavant_Values import dunavant_values as dunavant_values_import

import Dunavant_Values

class PolyDunavant():
    def __init__(self):
        self.order_values = self._initialize_order_values()
    
    @staticmethod
    def _initialize_order_values():
        dunavant_values = dunavant_values_import
        for order, values in dunavant_values.items():
            w, alpha, beta, gamma = PolyDunavant._apply_rotations(values)
            values['w'], values['alpha'], values['beta'], values['gamma'] = w, alpha, beta, gamma
        return dunavant_values
    
    @staticmethod
    def _apply_rotations(values):
        w, alpha, beta, gamma, degeneracy = (
            np.array(values['w']), 
            np.array(values['alpha']), 
            np.array(values['beta']), 
            np.array(values['gamma']), 
            values['degeneracy']
        )

        # Start with the initial values (accounting for the first set of values without rotation)
        rotated_weights = [w[0]]
        rotated_alpha = [alpha[0]]
        rotated_beta = [beta[0]]
        rotated_gamma = [gamma[0]]
        
        # Handle the first set of rotations based on the first level of degeneracy
        for i in range(degeneracy[1]):
            idx = i + degeneracy[0]  # Adjust for the index offset
            
            # Append the original rotation
            rotated_weights.append(w[idx])
            rotated_alpha.append(alpha[idx])
            rotated_beta.append(beta[idx])
            rotated_gamma.append(gamma[idx])

            # Append rotations
            rotations = [(alpha[idx], beta[idx], gamma[idx]),  # Original
                        (gamma[idx], alpha[idx], beta[idx]),  # Rotation 1
                        (beta[idx], gamma[idx], alpha[idx])]  # Rotation 2
            for a, b, g in rotations[1:]:  # Exclude the original since it's already added
                rotated_weights.append(w[idx])  # Weight is the same for each rotation
                rotated_alpha.append(a)
                rotated_beta.append(b)
                rotated_gamma.append(g)

        # Handle the second set of rotations based on the second level of degeneracy
        for i in range(degeneracy[2]):
            idx = i + degeneracy[1] + degeneracy[0]  # Adjust for the index offset
            
            # For the second level of degeneracy, all six permutations need to be considered
            permutations = [(alpha[idx], beta[idx], gamma[idx]), 
                            (gamma[idx], alpha[idx], beta[idx]), 
                            (beta[idx], gamma[idx], alpha[idx]),
                            (alpha[idx], gamma[idx], beta[idx]), 
                            (gamma[idx], beta[idx], alpha[idx]), 
                            (beta[idx], alpha[idx], gamma[idx])]
            for a, b, g in permutations:
                rotated_weights.append(w[idx])  # Weight is the same for each permutation
                rotated_alpha.append(a)
                rotated_beta.append(b)
                rotated_gamma.append(g)
        
        # Convert lists to numpy arrays before returning
        return (np.array(rotated_weights), np.array(rotated_alpha), 
                np.array(rotated_beta), np.array(rotated_gamma))
    def get_values(self, P):
        if P in self.order_values:
            order = self.order_values[P]
            return (order['w'], order['alpha'], order['beta'], order['gamma'], order['ng'])
        else:
            print("Order not precomputed. Please add the order values to _initialize_order_values.")
            return np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1), 0
        

# P=4
# optimized_poly_dun = PolyDunavant()
# weight2, alpha2, beta2, gamma2, ng2 = optimized_poly_dun.get_values(P=P)