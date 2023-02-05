import numpy as np
import meshing as msh

E = 1000 #Young's modulus in MPa
v = 0.25   #Poisson's ratio
pi = np.pi

R = 1  # Radius of circle
N = 8  # number of division in r direction
theta = np.pi/8  # number of division in theta direction
mps = 2*N-1  # number of mesh per anglular division 
sec = int(2*np.pi/theta)  # number of angular division
nodenum = N*sec + 1  # number of node
meshnum = mps*sec  # number of mesh

p = np.zeros([2, nodenum])  # Define matrix for coordinate
M = np.zeros([meshnum, 3])  # Define matrix for mesh

msh.Findcoor(p, R, N, theta)
msh.FindIndex(M, N, theta)
M = M.astype(int)
# msh.Mesh(p, M)

N0 = [0,0]
N1 = [1,0]
N2 = [1,2]
N3 = [0,2]

pp = [N0, N1, N2, N3]

M1 = [0, 1, 2]
M2 = [0, 2, 3]

MM = np.array([M1, M2]) 


class K_pstrain:
    def __init__(self, M, P,  E, v):
        self.E = E
        self.v = v
        self.P = np.transpose(P)
        self.x1 = self.P[M[0]][0]
        self.y1 = self.P[M[0]][1]
        self.x2 = self.P[M[1]][0]
        self.y2 = self.P[M[1]][1]
        self.x3 = self.P[M[2]][0]
        self.y3 = self.P[M[2]][1]
    
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

def Ktot(M, p, E, v):
    Ksize = 2*len(np.transpose(p))
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


def Kreduce(K, N, sec):
    rk = np.delete(K, int(2*N*(sec/2 + 1) + 1), 0)
    rk = np.delete(rk, int(2*N*(sec/2 + 1) + 1), 1)
    for i in range(int(2*N*(sec/2 + 1)), int(2*(1+N*sec/2)-2), -2):
        rk = np.delete(rk,i,0)
        rk = np.delete(rk,i,1)
    for i in range(2*N, -2, -2):
        rk = np.delete(rk,i,0)
        rk = np.delete(rk,i,1)
    return rk

F = np.zeros([2*nodenum,1])
F[8][0] = 0
F[9][0] = -500
def Freduce(F, N, sec):
    fk = np.delete(F, int(2*N*(sec/2 + 1) + 1), 0)
    for i in range(int(2*N*(sec/2 + 1)), int(2*(1+N*sec/2)-2), -2):
        fk = np.delete(fk,i,0)
    for i in range(2*N, -2, -2):
        fk = np.delete(fk,i,0)
    return fk

def cal_U(K, F):
    invK = np.linalg.inv(K)
    u = invK.dot(F)
    
    return u

def Uremake(u, N, sec):
    z = np.zeros([1,1])
    U = np.insert(u,0,z,0)
    for i in range(2, 2*N +2, 2):
        U = np.insert(U,i,z,0)
    for i in range(int(2*(1+N*sec/2)), int(2*N*(sec/2 + 1) + 2),  2):
        U = np.insert(U,i,z,0)
    U = np.insert(U, int(2*N*(sec/2 + 1) + 1), z, 0)
    
    return U

def Deform(p, u):
    p = np.transpose(p)
    for i in range(len(p)):
        p[i][0] += u[2*i][0]

        if p[i][1] + u[2*i + 1][0] < -1:
            p[i][1] = -1
        else:
            p[i][1] += u[2*i + 1][0]

    return np.transpose(p)

# K1 = K_pstress(M1, p, E, v).K()
# K2 = K_pstrain(M2, p, E, v).K()
TK = Ktot(M, p , E, v )
RK = Kreduce(TK, N, sec)
RF = Freduce(F, N, sec)
u = cal_U(RK, RF)
uu = Uremake(u, N, sec)
pp = Deform(p,uu)
# RK = np.delete(TK, 0, 0)
print(len(np.transpose(p)))
print(np.shape(TK))
print(np.shape(RK))
# print(np.shape(uu))
# print(uu)

msh.Mesh(pp, M)
# print(len(TK))
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

