import numpy as np

def E_field_essentials(k, direction, polarization_angle):
    # Calculate wave vector components
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

    # Cross products
    Vect1 = np.cross([kx_phi, ky_phi, kz_phi], k_vector)
    Vect2 = np.cross([kx_theta, ky_theta, kz_theta], k_vector)

    # Normalize the perpendicular vectors
    pol_vect1 = Vect1 / np.linalg.norm(Vect1)
    pol_vect2 = Vect2 / np.linalg.norm(Vect2)

    polarization_angle = np.pi/4
    final_polarization = np.sin(polarization_angle)*pol_vect2 + np.cos(polarization_angle)*pol_vect1
    final_polarization = final_polarization / np.linalg.norm(final_polarization)
    return k_vector, final_polarization

def E_field_in_point(pos, k_vector, polarization, amp):
    phase_shift = np.exp(-1j*np.dot(k_vector,pos))
    E_field = amp * polarization * phase_shift
    return E_field