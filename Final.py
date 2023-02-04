import numpy as np
import scipy as sp
import sys


E = 1000 #Young's modulus in MPa
v = 0.25   #Poisson's ratio

N0 = [0,0]
N1 = [1,0]
N2 = [1,2]
N3 = [0,2]

p = [N0, N1, N2, N3]

M1 = [0, 1, 2]
M2 = [0, 2, 3]

M = np.array([M1, M2]) 

class K_pstrain:
    def __init__(self, M, P,  E, v):
        self.E = E
        self.v = v
        self.P = P
        self.x1 = P[M[0]][0]
        self.y1 = P[M[0]][1]
        self.x2 = P[M[1]][0]
        self.y2 = P[M[1]][1]
        self.x3 = P[M[2]][0]
        self.y3 = P[M[2]][1]
    
    def B(self):
        b = np.array([[self.y2 - self.y3,           0        ,   self.y3 - self.y1,           0        ,   self.y1 - self.y2,           0        ], \
                      [        0        ,   self.x3 - self.x2,           0        ,   self.x1 - self.x3,           0        ,   self.x2 - self.x1], \
                      [self.x3 - self.x2,   self.y2 - self.y3,   self.x1 - self.x3,   self.y3 - self.y1,   self.x2 - self.x1,   self.y1 - self.y2]])
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
    def __init__(self, M, P, E, v):
        super().__init__(M, P, E, v)
    
    def D(self):
        d = np.array([[   1   , self.v,      0      ], \
                      [ self.v,   1   ,      0      ], \
                      [   0   ,   0   , (1-self.v)/2]])

        D = E / (1+v) / (1-v) * d

        return D


def totK(M, p, E, v):
    Ksize = 2*len(p)
    K = np.zeros([Ksize, Ksize])
    
    # Making stiffness matrix
    for i in range(len(M)):
        Ki = K_pstress(M[i], p, E, v).K()

        for j in range(len(M[i])):
            idr = M[i][j]
            for k in range(len(M[i])):
                idc = M[i][k]
                
                K[2*idr][2*idc] += Ki[2*j][2*k]
                K[2*idr][2*idc+1] += Ki[2*j][2*k+1]
                K[2*idr+1][2*idc] += Ki[2*j+1][2*k]
                K[2*idr+1][2*idc+1] += Ki[2*j+1][2*k+1]

    return K



    # for row in range(len(Ki)):
    #     for col in range(len(Ki)):
    #         K[col][row] = Ki[M[0]]
    
    
    return K

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


def reduce(K):
    rk = np.zeros([4,4])
    rk[0] = [K[2][2], K[2][4], K[2][5], K[2][7]]
    rk[1] = [K[4][2], K[4][4], K[4][5], K[4][7]]
    rk[2] = [K[5][2], K[5][4], K[5][5], K[5][7]]
    rk[3] = [K[7][2], K[7][4], K[7][5], K[7][7]] 
    return rk

K1 = K_pstress(M1, p, E, v).K()
K2 = K_pstrain(M2, p, E, v).K()
TK = totK(M, p , E, v )
print(np.around(TK,2))
# KK = totK(K1,K2)
# RK = reduce(KK)
# iK = np.linalg.inv(KK)
# U = cal_U(RK, f)

# F1 = cal_F(reduce(K1), U)
# F2 = cal_F(reduce(K2), U)
# print(np.around(KK,2))
# print(np.around(RK,2))
# print(U)
# # print(np.around(U,2))
# print(np.around(KK.dot(iK),2))

