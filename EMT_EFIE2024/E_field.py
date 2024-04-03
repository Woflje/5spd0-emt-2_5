
import numpy as np

def E_field_essentials(k, direction, polarization_angle):
    # Calculate wave vector components, and wave vector
    kx = k * np.sin(direction[0]) * np.cos(direction[1])
    ky = k * np.sin(direction[0]) * np.sin(direction[1])
    kz = k * np.cos(direction[0])
    k_vector = np.array([kx, ky, kz])

    # Slightly altered vectors for phi and theta
    kx_phi = k * np.sin(direction[0]+0.001) * np.cos(direction[1])
    ky_phi = k * np.sin(direction[0]+0.001) * np.sin(direction[1])
    kz_phi = k * np.cos(direction[0]+0.001)

    kx_theta = k * np.sin(direction[0]) * np.cos(direction[1]+0.004)
    ky_theta = k * np.sin(direction[0]) * np.sin(direction[1]+0.004)
    kz_theta = k * np.cos(direction[0])

    Vect1 = np.cross([kx_phi, ky_phi, kz_phi], k_vector)
    Vect2 = np.cross([kx_theta, ky_theta, kz_theta], k_vector)

    pol_vect1 = Vect1 / np.linalg.norm(Vect1)
    pol_vect2 = Vect2 / np.linalg.norm(Vect2)

    final_polarization = np.sin(polarization_angle)*pol_vect2 + np.cos(polarization_angle)*pol_vect1
    final_polarization = final_polarization / np.linalg.norm(final_polarization)
    return k_vector, final_polarization

def E_field(vertices, rwgs, areas, dunavant_positions, k, amp, direction, polarization_angle, dunavant_weight):
    [k_vector, polarization] = E_field_essentials(k,direction,polarization_angle)

    edge_vectors = np.subtract(vertices[rwgs[:,0]], vertices[rwgs[:,1]])
    length_common_edges = np.linalg.norm(edge_vectors, axis=1)

    # Calculate areas and select traingle 1 and 2 vertices for all RWGs
    A1 = areas[rwgs[:,4]]
    A2 = areas[rwgs[:,5]]
    T1 = vertices[rwgs[:,0:3]]
    T2 = vertices[rwgs[:,[0,1,3]]]

    # Dunavant positions for the two triangles in all RWGs
    dunavant_positions1 = dunavant_positions[:,0]
    dunavant_positions2 = dunavant_positions[:,1]

    dot_product1 = np.einsum('i,npi->np', k_vector, dunavant_positions1)
    dot_product2 = np.einsum('i,npi->np', k_vector, dunavant_positions2)

    temp1 = amp * polarization[:, np.newaxis, np.newaxis] * np.exp(-1j * dot_product1)
    temp2 = amp * polarization[:, np.newaxis, np.newaxis] * np.exp(-1j * dot_product2)

    diff_positions1 = dunavant_positions1 - T1[:,2][:, np.newaxis, :]
    diff_positions2 = dunavant_positions2 - T2[:,2][:, np.newaxis, :]
    scaling_factor1 = length_common_edges[:, np.newaxis] / (2 * A1[:, np.newaxis])
    scaling_factor2 = length_common_edges[:, np.newaxis] / (2 * A2[:, np.newaxis])
    temp1 = (amp * np.exp(-1j * dot_product1))[..., np.newaxis] * polarization
    temp2 = (amp * np.exp(-1j * dot_product2))[..., np.newaxis] * polarization

    temp2_1 = np.einsum('npi,npi->np', diff_positions1 * scaling_factor1[:, :, np.newaxis], temp1)
    temp2_2 = np.einsum('npi,npi->np', diff_positions2 * scaling_factor2[:, :, np.newaxis], temp2)

    E_t1 = A1 * np.sum(temp2_1 * dunavant_weight[np.newaxis,:], axis=1)
    E_t2 = A2 * np.sum(temp2_2 * dunavant_weight[np.newaxis,:], axis=1)

    E = E_t1 - E_t2

    return E.reshape(-1, 1)