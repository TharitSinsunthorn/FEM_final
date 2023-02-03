import numpy as np
import sys


E = 1000 #Young's modulus in MPa
v = 0.25   #Poisson's ratio

N1 = [0,0]
N2 = [1,0]
N3 = [1,1]

M = [N1, N2, N3]

def mkB(M):
    B = np.array([[M[1][1]-M[2][1],        0       , M[2][1]-M[0][1],        0       , M[0][1]-M[1][1],        0       ], \
                  [       0       , M[2][0]-M[1][0],        0       , M[0][0]-M[2][0],        0       , M[1][0]-M[0][0]], \
                  [M[2][0]-M[1][0], M[1][1]-M[2][1], M[0][0]-M[2][0], M[2][1]-M[0][1], M[1][0]-M[0][0], M[0][1]-M[1][1]]])
    return B

def mkD_pstrain(E, v):
    d = np.array([[1-v,  v ,     0    ], \
                 [ v , 1-v,     0    ], \
                 [ 0 ,  0 , (1-2*v)/2]])

    D = E / (1+v) / (1-2*v) * d

    return D




B = mkB(M)
D = mkD_pstrain(E,v)

print(B)
print(D)
