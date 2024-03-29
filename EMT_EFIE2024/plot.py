import numpy as np
import matplotlib.pyplot as plt

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
