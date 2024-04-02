import numpy as np
from integrators.dunavant.dunavant_values import dunavant_values

def apply_rotations_new(values):
    w, alpha, beta, gamma, degeneracy = (
        np.array(values['w']), 
        np.array(values['alpha']), 
        np.array(values['beta']), 
        np.array(values['gamma']), 
        values['degeneracy']
    )

    rotated_alpha = []
    rotated_beta = []
    rotated_gamma = []
    rotated_weights = []
    
    rotated_alpha.append(alpha[0])
    rotated_beta.append(beta[0])
    rotated_gamma.append(gamma[0])
    rotated_weights.append(w[0])
    
    for y in range(degeneracy[1]):
        it = y+degeneracy[0]
        rotated_alpha.append(alpha[it])
        rotated_beta.append(beta[it])
        rotated_gamma.append(gamma[it])
        rotated_weights.append(w[it])
        
        rotated_alpha.append(gamma[it])
        rotated_beta.append(alpha[it])
        rotated_gamma.append(beta[it])
        rotated_weights.append(w[it])
        
        rotated_alpha.append(beta[it])
        rotated_beta.append(gamma[it])
        rotated_gamma.append(alpha[it])
        rotated_weights.append(w[it])
        
    for z in range(degeneracy[2]):
        it = z+degeneracy[1]+degeneracy[0]
        rotated_alpha.append(alpha[it])
        rotated_beta.append(beta[it])
        rotated_gamma.append(gamma[it])
        rotated_weights.append(w[it])
        
        rotated_alpha.append(gamma[it])
        rotated_beta.append(alpha[it])
        rotated_gamma.append(beta[it])
        rotated_weights.append(w[it])
        
        rotated_alpha.append(beta[it])
        rotated_beta.append(gamma[it])
        rotated_gamma.append(alpha[it])
        rotated_weights.append(w[it])
        
        rotated_alpha.append(alpha[it])
        rotated_beta.append(gamma[it])
        rotated_gamma.append(beta[it])
        rotated_weights.append(w[it])
        
        rotated_alpha.append(gamma[it])
        rotated_beta.append(beta[it])
        rotated_gamma.append(alpha[it])
        rotated_weights.append(w[it])
        
        rotated_alpha.append(beta[it])
        rotated_beta.append(alpha[it])
        rotated_gamma.append(gamma[it])
        rotated_weights.append(w[it])
    
    return (np.array(rotated_weights), np.array(rotated_alpha), 
            np.array(rotated_beta), np.array(rotated_gamma))

def apply_rotations(values):
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

def initialize_order_values():
    for order, values in dunavant_values.items():
        w, alpha, beta, gamma = apply_rotations_new(values)
        values['w'], values['alpha'], values['beta'], values['gamma'] = w, alpha, beta, gamma
    return dunavant_values

def get_dunavant_values(order):
    values = dunavant_values[order]
    return np.vstack([values['alpha'],values['beta'],values['gamma']]).T, values['w']

dunavant_values = initialize_order_values()