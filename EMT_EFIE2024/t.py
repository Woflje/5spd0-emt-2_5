import numpy as np

# Assuming the required variables and arrays (vertices, rwgs, dunavant_positions, areas, k_vector, amp, polarization, dunavant_weight) are predefined

# Calculate the lengths of common edges for all RWGs
edge_vectors = np.subtract(vertices[rwgs[:,0]], vertices[rwgs[:,1]])
length_common_edges = np.linalg.norm(edge_vectors, axis=1)

# Calculate areas and select T1 and T2 vertices for all RWGs
A1 = areas[rwgs[:,4]]
A2 = areas[rwgs[:,5]]
T1 = vertices[rwgs[:,0:3]]
T2 = vertices[rwgs[:,[0,1,3]]]

# Dunavant positions for the two triangles in all RWGs
dunavant_positions1 = dunavant_positions[:,0]
dunavant_positions2 = dunavant_positions[:,1]

# Dot product for all positions in dunavant_positions1 and dunavant_positions2 with k_vector
dot_product1 = np.einsum('i,npi->np', k_vector, dunavant_positions1)
dot_product2 = np.einsum('i,npi->np', k_vector, dunavant_positions2)

# Compute 'temp' for all positions at once for both sets of Dunavant positions
temp1 = amp * polarization[:, np.newaxis, np.newaxis] * np.exp(-1j * dot_product1)
temp2 = amp * polarization[:, np.newaxis, np.newaxis] * np.exp(-1j * dot_product2)

# Compute 'temp2' for all positions at once for both sets of Dunavant positions
diff_positions1 = dunavant_positions1 - T1[:,2][:, np.newaxis, :]
diff_positions2 = dunavant_positions2 - T2[:,2][:, np.newaxis, :]
scaling_factor1 = length_common_edges[:, np.newaxis] / (2 * A1[:, np.newaxis])
scaling_factor2 = length_common_edges[:, np.newaxis] / (2 * A2[:, np.newaxis])
temp1 = amp * np.exp(-1j * dot_product1)  # This should now be of shape (N, P), matching the dot_product1 shape
temp2 = amp * np.exp(-1j * dot_product2)  # Similarly, for dot_product2

# Adjust polarization application to each temp
temp1 = temp1[..., np.newaxis] * polarization  # Broadcasting polarization, shape becomes (N, P, 3)
temp2 = temp2[..., np.newaxis] * polarization  # Broadcasting polarization, shape becomes (N, P, 3)

# Now, we ensure diff_positions and temps are aligned for einsum operation
# The operation on diff_positions*scaling_factor does not change, as it already produces an array shaped (N, P, 3)
temp2_1 = np.einsum('npi,npi->np', diff_positions1 * scaling_factor1[:, :, np.newaxis], temp1)
temp2_2 = np.einsum('npi,npi->np', diff_positions2 * scaling_factor2[:, :, np.newaxis], temp2)

# Compute 'E_t1' and 'E_t2' by summing over all 'temp2' elements, scaled by 'A1', 'A2' and 'dunavant_weight'
E_t1 = A1 * np.sum(temp2_1 * dunavant_weight[np.newaxis,:], axis=1)  # Ensure dunavant_weight is broadcasted correctly
E_t2 = A2 * np.sum(temp2_2 * dunavant_weight[np.newaxis,:], axis=1)

# Now E_t1 and E_t2 should be of shape (N,)
# Calculate E2 as the difference between E_t1 and E_t2
E2 = E_t1 - E_t2

# Ensure E2 is shaped (N,), or reshape if you specifically need (N,1)
E2 = E2.reshape(-1, 1)  # Reshapes E2 to (N, 1) if required
