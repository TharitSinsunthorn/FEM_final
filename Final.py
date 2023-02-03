import numpy as np
import scipy as sp
import sys


E = 1000 #Young's modulus in MPa
v = 0.25   #Poisson's ratio

N1 = [0,0]
N2 = [1,0]
N3 = [1,2]
N4 = [0,2]

M1 = [N1, N2, N3]
M2 = [N1, N3, N4]

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
        b = np.array([[self.M[1][1]-self.M[2][1],             0            , self.M[2][1]-self.M[0][1],             0            , self.M[0][1]-self.M[1][1],             0            ], \
                      [            0            , self.M[2][0]-self.M[1][0],             0            , self.M[0][0]-self.M[2][0],             0            , self.M[1][0]-self.M[0][0]], \
                      [self.M[2][0]-self.M[1][0], self.M[1][1]-self.M[2][1], self.M[0][0]-self.M[2][0], self.M[2][1]-self.M[0][1], self.M[1][0]-self.M[0][0], self.M[0][1]-self.M[1][1]]])
        B = 1/(2*self.mesharea()) * b

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

class K_pstress(K_pstrain):
    def __init__(self, M, E, v):
        super().__init__(M, E, v)
    
    def D(self):
        d = np.array([[   1   , self.v,      0      ], \
                      [ self.v,   1   ,      0      ], \
                      [   0   ,   0   , (1-self.v)/2]])

        D = E / (1+v) / (1-v) * d

        return D


u = np.zeros([6,1])
u[4] = 0.01
u[5] = 0

def cal_F(K, U):
    F = K.dot(U)

    return F


f = np.array([[0],[0],[1],[1]])
def cal_U(K, F):
    invK = np.linalg.inv(K)
    U = invK.dot(F)
    
    return U

def combineK(K1, K2):
    tk = np.zeros([8,8])
    tk[0] = [K1[0][0]+K2[0][0], K1[0][1]+K2[0][1], K1[0][2], K1[0][3], K1[0][4]+K2[0][2], K1[0][5]+K2[0][3], K2[0][4], K2[0][5]] 
    tk[1] = [K1[1][0]+K2[1][0], K1[1][1]+K2[1][1], K1[1][2], K1[1][3], K1[1][4]+K2[1][2], K1[1][5]+K2[1][3], K2[1][4], K2[1][5]] 
    tk[2] = [K1[2][0]         , K1[2][1]         , K1[2][2], K1[2][3], K1[2][4]         , K1[2][5]         ,     0   ,    0    ] 
    tk[3] = [K1[3][0]         , K1[3][1]         , K1[3][2], K1[3][3], K1[3][4]         , K1[3][5]         ,     0   ,    0    ] 
    tk[4] = [K1[4][0]+K2[2][0], K1[4][1]+K2[2][1], K1[4][2], K1[4][3], K1[4][4]+K2[2][2], K1[4][5]+K2[2][3], K2[2][4], K2[2][5]] 
    tk[5] = [K1[5][0]+K2[3][0], K1[5][1]+K2[3][1], K1[5][2], K1[5][3], K1[5][4]+K2[3][2], K1[5][5]+K2[3][3], K2[3][4], K2[3][5]] 
    tk[6] = [         K2[4][0],          K2[4][1],     0   ,     0   ,          K2[4][2],          K2[4][3], K2[4][4], K2[4][5]] 
    tk[7] = [         K2[5][0],          K2[5][1],     0   ,     0   ,          K2[5][2],          K2[5][3], K2[5][4], K2[5][5]] 
    

    return tk

def reduce(K):
    rk = np.zeros([4,4])
    rk[0] = [K[2][2], K[2][4], K[2][5], K[2][7]]
    rk[1] = [K[4][2], K[4][4], K[4][5], K[4][7]]
    rk[2] = [K[5][2], K[5][4], K[5][5], K[5][7]]
    rk[3] = [K[7][2], K[7][4], K[7][5], K[7][7]] 
    return rk

K1 = K_pstress(M1, E, v).K()
K2 = K_pstress(M2, E, v).K()

KK = combineK(K1,K2)
RK = reduce(KK)
iK = np.linalg.inv(KK)
U = cal_U(RK, f)
print(np.around(KK,2))
print(np.around(RK,2))
print(U)
# print(np.around(U,2))
# print(np.around(KK.dot(iK),2))

