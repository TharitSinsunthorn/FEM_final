import numpy as np
import sys


E = 1000 #Young's modulus in MPa
v = 0.25   #Poisson's ratio

N1 = [0,0]
N2 = [1,0]
N3 = [1,1]

M = [N1, N2, N3]

class K_pstrain:
    def __init__(self, M, E, v):
        self.E = E
        self.v = v
        self.M = M
        self.x1 = M[0][0]
        self.y1 = M[0][1]
        self.x2 = M[1][0]
        self.y2 = M[1][1]
        self.x3 = M[2][0]
        self.y3 = M[2][1]
    
    def B(self):
        B = np.array([[self.M[1][1]-self.M[2][1],             0            , self.M[2][1]-self.M[0][1],             0            , self.M[0][1]-self.M[1][1],             0            ], \
                      [            0            , self.M[2][0]-self.M[1][0],             0            , self.M[0][0]-self.M[2][0],             0            , self.M[1][0]-self.M[0][0]], \
                      [self.M[2][0]-self.M[1][0], self.M[1][1]-self.M[2][1], self.M[0][0]-self.M[2][0], self.M[2][1]-self.M[0][1], self.M[1][0]-self.M[0][0], self.M[0][1]-self.M[1][1]]])
        # B = np.array([[M[1][1]-M[2][1],        0       , M[2][1]-M[0][1],        0       , M[0][1]-M[1][1],        0       ], \
        #             [       0       , M[2][0]-M[1][0],        0       , M[0][0]-M[2][0],        0       , M[1][0]-M[0][0]], \
        #             [M[2][0]-M[1][0], M[1][1]-M[2][1], M[0][0]-M[2][0], M[2][1]-M[0][1], M[1][0]-M[0][0], M[0][1]-M[1][1]]])
        return B

    def D(self):
        d = np.array([[1-self.v,    self.v,        0      ], \
                      [ self.v ,  1-self.v,        0      ], \
                      [ 0      ,      0   , (1-2*self.v)/2]])

        D = E / (1+v) / (1-2*v) * d

        return D

    def mesharea(self):
        a = 0.5 * np.linalg.det(np.array([[1,self.x1,self.y1], [1,self.x2,self.y2], [1,self.x3,self.y3]]))
        return a

    def K(self):
        a = self.mesharea()
        B = self.B()
        D = self.D()
        Bt = np.transpose(B)

        k = a * Bt.dot(D.dot(B))
        
        return k

K = K_pstrain(M,E,v)

print(K.K())

