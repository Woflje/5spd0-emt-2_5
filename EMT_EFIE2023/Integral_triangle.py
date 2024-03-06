import numpy as np
import scipy


def area_triangle(T):
    l1 = np.linalg.norm(np.subtract(T[1],T[0]))
    l2 = np.linalg.norm(np.subtract(T[0],T[2]))
    l3 = np.linalg.norm(np.subtract(T[2],T[1]))
    
    A = 1/4*np.sqrt(4*l1**2*l2**2-(l1**2+l2**2-l3**2)**2)
    return A

def int_bastest(G,T1,T2):
    Int1 = 1
    
    return G(T1[1],T2[2])



def int_test(E,T):    
    l_edge = np.linalg.norm(np.subtract(T[0],T[1]))
    A = area_triangle(T)
    scalar = lambda r: np.dot(l_edge/(2*A)*np.subtract(r,T[2]),E(r))
    
    Itest = int_triangle(scalar,T)
    
    return Itest


def int_triangle(Phi,T):
    ng = 4;
    A_tot = area_triangle(T)
    
    temp = 0;
    weight = [-0.5625,0.520833333333333,0.520833333333333,0.520833333333333]
    alpha = [0.333333333333,0.6,0.2,0.2]
    beta = [0.333333333333,0.2,0.6,0.2]
    gamma = [0.333333333333,0.2,0.2,0.6]
    
    
    for i in range(ng):
        r = np.multiply(alpha[i],T[0]) + np.multiply(beta[i],T[1]) + np.multiply(gamma[i],T[2])
        temp = temp + Phi(r)*weight[i]
        
    return A_tot*temp

T1 = [[1,1,1],[1,5,3],[10,3,2]]
T2 = [[1,1,1],[1,5,3],[8,3,5]]
Phi = lambda r: np.linalg.norm(r)**2
E = lambda r: r
G_ft2 = lambda r,rp: [np.add(r,rp),np.add(r,rp),np.add(r,rp)]

A = int_triangle(Phi,T1)

r1 = int_test(E,T1)

val = int_bastest(G, T1, T2)
print(val)

    
    