import numpy as np
from FastNP import fastNorm

def int_triangle(integrand,A,dunavant_pos,dunavant_weight):
    phi_t=np.array([integrand(r) for r in dunavant_pos])
    return A*sum(np.multiply(phi_t[:,np.newaxis],dunavant_weight[:,np.newaxis]))[0]

    #Function that calculates the integral over a triangle for a vector field
def int_test(integrand,triangle_points,A,dunavant_pos,dunavant_weight):
    l_edge = fastNorm(np.subtract(triangle_points[0],triangle_points[1])) #Length of common edge
    scalar = lambda r: np.dot(l_edge/(2*A)*np.subtract(r,triangle_points[2]),integrand(r)) #This function ensures that we only have to integrate over a scalar function as the test/basis function dotted with the vector field results in a scalar function.
    Itest = int_triangle(scalar,A,dunavant_pos,dunavant_weight) #Integrate the new scalar function using the int_triangle function.
    return Itest