import numpy as np
import matplotlib.pyplot as plt
import FastNP as fnp

def E_in_vectorized2(amp, dir, pos, polarization, k):
    # Convert spherical to cartesian coordinates for the wave direction
    kx = k * np.sin(dir[0]) * np.cos(dir[1])
    ky = k * np.sin(dir[0]) * np.sin(dir[1])
    kz = k * np.cos(dir[0])
    k_vec = np.array([kx, ky, kz])

    # Calculate electric field vector
    # Assuming polarization affects the amplitude direction, simplified here
    E0 = amp * np.array([np.cos(polarization), np.sin(polarization), 0])  # Simplified polarization
    
    # Calculate phase shift for each position
    phase_shift = np.exp(-1j * np.dot(pos, k_vec))

    # Apply amplitude and polarization
    E_field = E0[:, np.newaxis] * phase_shift
    
    return E_field

def E_in_vectorized(amp, dir, pos, polarization, wavelength):
    # Determine wavenumber
    k = 2*np.pi / wavelength

    #Determine wave vector using input angles phi (dir[0]) and theta (dir[1])
    kx = k*np.sin(dir[0])*np.cos(dir[1])
    ky = k*np.sin(dir[0])*np.sin(dir[1])
    kz = k*np.cos(dir[0])
    wave_vector = np.array([kx, ky, kz])

    #Find values for a vector with a slightly different phi value
    kx_phi = k*np.sin(dir[0]+0.001)*np.cos(dir[1])
    ky_phi = k*np.sin(dir[0]+0.001)*np.sin(dir[1])
    kz_phi = k*np.cos(dir[0]+0.001)
    phi_vect = np.array([kx_phi, ky_phi, kz_phi])

    #Find values for a vector with a slightly different theta value
    kx_th = k*np.sin(dir[0])*np.cos(dir[1]+0.004)
    ky_th = k*np.sin(dir[0])*np.sin(dir[1]+0.004)
    kz_th = k*np.cos(dir[0])
    theta_vect = np.array([kx_th, ky_th, kz_th])

    # Cross product of the wave vector and the changed phi and theta vectors to find perpendicular vectors
    Vect1 = fnp.fastCross(phi_vect, wave_vector)
    Vect2 = fnp.fastCross(theta_vect, wave_vector)

    # Normalize the perpendicular vectors
    pol_vect1 = Vect1 / fnp.fastNorm(Vect1)
    pol_vect2 = Vect2 / fnp.fastNorm(Vect2)

    # Allow user to choose polarization by using input angle between the two normalized polarization vectors
    polarization_vect = np.sin(polarization) * pol_vect2 + np.cos(polarization) * pol_vect1
    polarization_vect = polarization_vect / fnp.fastNorm(polarization_vect)

    # Dot product to check if polarization is perpendicular to the wave vector
    if np.dot(wave_vector, polarization_vect) > 1e-10:
        print("The polarization and wavenumber are not perpendicular!")

    # Calculate the E-field for each position
    phase_shift = np.exp(-1j * np.dot(wave_vector, pos.T))  # Dot product and phase shift for each position
    E_field = amp * polarization_vect[:, np.newaxis] * phase_shift  # Apply amplitude and polarization

    return E_field.T  # Return transpose to match the expected shape (N, 3)

def E_in(amp, dir, pos, polarization, wavelength): # Function to calculate the incoming electric field
    #Determine wavenumber
    k = 2 * np.pi / (wavelength)
    #Determine wave vector using input angles phi (dir[0]) and theta (dir[1])
    kx = k*np.sin(dir[0])*np.cos(dir[1])
    ky = k*np.sin(dir[0])*np.sin(dir[1])
    kz = k*np.cos(dir[0])

    #Find values for a vector with a slightly different phi value
    kx_phi = k*np.sin(dir[0]+0.001)*np.cos(dir[1])
    ky_phi = k*np.sin(dir[0]+0.001)*np.sin(dir[1])
    kz_phi = k*np.cos(dir[0]+0.001)
    #Cross product of the wave vector and the changed phi vector to find a perpendicular vector
    Vect1 = fnp.fastCross(np.array([kx_phi,ky_phi,kz_phi]), np.array([kx,ky,kz]))

    #Find values for a vector with a slightly different theta value
    kx_th = k*np.sin(dir[0])*np.cos(dir[1]+0.004)
    ky_th = k*np.sin(dir[0])*np.sin(dir[1]+0.004)
    kz_th = k*np.cos(dir[0])
    #Cross product of the wave vector and the changed theta vector to find a perpendicular vector
    Vect2 = fnp.fastCross(np.array([kx_th,ky_th,kz_th]), np.array([kx,ky,kz])) 

    #Normalize the perpendicular vectors
    pol_vect1 = np.multiply((1./fnp.fastNorm(Vect1)),Vect1)
    #print('k phi')
    #print(pol_vect1)
    pol_vect2 = np.multiply((1./fnp.fastNorm(Vect2)),Vect2)
    #print('k theta')
    #print(pol_vect2)
    #Allow user to choose polarization by using input angle between the two normalized polarization vectors
    polarization = np.sin(polarization)*pol_vect2+np.cos(polarization)*pol_vect1
    polarization = np.multiply((1/fnp.fastNorm(polarization)),polarization)
    #print('polarization')
    #print(polarization)

    #Check if polarization is perpendicular to the wave vector
    if np.dot([kx, ky, kz], polarization) > 1e-10:
        print("The polarization and wavenumber are not perpendicular!")

    #return the full plane wave
    #return np.array([1,0,0])
    return np.multiply(amp*polarization,np.exp(np.multiply(-1J, np.dot([kx, ky, kz],pos))))


def check_sing(vertice11, vertice12, vertice13, vertice21, vertice22, vertice23):  # Function to check if the two triangles for integration are singular
    # Check each vertice if there is an overlapping vertice of the second triangle
    if np.array_equal(vertice11,vertice21) or np.array_equal(vertice11, vertice22) or np.array_equal(vertice11, vertice23):
        # The triangles are Singular
        return True
    elif np.array_equal(vertice12,vertice21) or np.array_equal(vertice12, vertice22) or np.array_equal(vertice12, vertice23):
        # The triangles are Singular
        return True
    elif np.array_equal(vertice13,vertice21) or np.array_equal(vertice13, vertice22) or np.array_equal(vertice13, vertice23):
        # The triangles are Singular
        return True
    else:
        # The triangles are not Singular
        return False

def check_sing_vectorized(vertice1, vertice2):
    for v1 in vertice1:
        for v2 in vertice2:
            if np.array_equal(vertice1, vertice2):
                return True
    return False

def plot_phi(rad,E_farfield):  # Function to plot the farfield plot for a constant theta
    # Create polar plot to plot in
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='polar')

    # Plot all farfield values for each angle rad of phi
    ax.plot(rad,np.real(E_farfield[0,:,0,0]))

    # Set labels and rotate the plot to match CST
    ax.set_xlabel('phi')
    ax.set_theta_zero_location("N")
    plt.show()

    return


def plot_theta(rad,E_farfield):  # Function to plot the farfield plot for a constant phi
    # Create polar plot to plot in
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='polar')

    # Plot all farfield values for each angle rad of theta
    plt.plot(rad,np.real(E_farfield[0,0,:,0]))

    # Set labels and rotate the plot to match CST
    ax.set_xlabel('theta')
    ax.set_theta_zero_location("N")
    plt.show()

    return
