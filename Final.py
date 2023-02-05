import numpy as np
import meshing as msh

E = 10000  #Young's modulus in MPa
v = 0.25   #Poisson's ratio
pi = np.pi

R = 1                           # Radius of circle
N = 4                           # number of division in r direction
theta = np.pi/4                 # number of division in theta direction
mps = 2*N-1                     # number of mesh per anglular division 
sec = int(2*np.pi/theta)        # number of angular division
nodenum = N*sec + 1             # number of node
meshnum = mps*sec               # number of mesh

p = np.zeros([2, nodenum])      # Define matrix for coordinate
M = np.zeros([meshnum, 3])      # Define matrix for mesh


msh.Findcoor(p, R, N, theta)    ## Find coordinate of each node
msh.FindIndex(M, N, theta)      ## Define index for all nodes to each mesh
M = M.astype(int)
# msh.Mesh(p, M)

## Apply the force matrix
F = np.zeros([2*nodenum,1])
F[2*N+1][0] = -500

## Class for calculate K matrix with plane strain condition
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
    
    # Calculate B matrix
    def B(self):
        b = np.array([[self.y2 - self.y3,           0        ,   self.y3 - self.y1,           0        ,   self.y1 - self.y2,           0        ], \
                      [        0        ,   self.x3 - self.x2,           0        ,   self.x1 - self.x3,           0        ,   self.x2 - self.x1], \
                      [self.x3 - self.x2,   self.y2 - self.y3,   self.x1 - self.x3,   self.y3 - self.y1,   self.x2 - self.x1,   self.y1 - self.y2]])
        B = 1/(2*self.mesharea()) * b

        return B

    # Calculate D matrix
    def D(self):
        d = np.array([[1-self.v,    self.v,        0      ], \
                      [ self.v ,  1-self.v,        0      ], \
                      [ 0      ,      0   , (1-2*self.v)/2]])

        D = E / (1+v) / (1-2*v) * d

        return D

    # Calculate area of mesh
    def mesharea(self):
        a = 0.5 * np.linalg.det(np.array([[1,self.x1,self.y1], [1,self.x2,self.y2], [1,self.x3,self.y3]]))
        return a

    # Calculate K matrix
    def K(self):
        a = self.mesharea()
        B = self.B()
        D = self.D()
        Bt = np.transpose(B)

        k = a * Bt.dot(D.dot(B))
        
        return k

## Class for calculate K matrix with plane stress condition
class K_pstress(K_pstrain):
    def __init__(self, M, P, E, v):
        super().__init__(M, P, E, v)
    
    # Calculate D matrix
    def D(self):
        d = np.array([[   1   , self.v,      0      ], \
                      [ self.v,   1   ,      0      ], \
                      [   0   ,   0   , (1-self.v)/2]])

        D = E / (1+v) / (1-v) * d

        return D

## Compute total stiffness of plane strain condition
def Ktot_pstrain(M, p, E, v):
    Ksize = 2*len(np.transpose(p))
    K = np.zeros([Ksize, Ksize])
    
    # Making stiffness matrix
    for i in range(len(M)):
        Ki = K_pstrain(M[i], p, E, v).K()

        for j in range(len(M[i])):
            idr = M[i][j]
            for k in range(len(M[i])):
                idc = M[i][k]
                
                K[2*idr][2*idc] += Ki[2*j][2*k]
                K[2*idr][2*idc+1] += Ki[2*j][2*k+1]
                K[2*idr+1][2*idc] += Ki[2*j+1][2*k]
                K[2*idr+1][2*idc+1] += Ki[2*j+1][2*k+1]

    return K

## Compute total stiffness of plane stress condition
def Ktot_pstress(M, p, E, v):
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
 
## Reduce dimension of K matrix
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

## Reduce dimension of F matrix
def Freduce(F, N, sec):
    fk = np.delete(F, int(2*N*(sec/2 + 1) + 1), 0)
    for i in range(int(2*N*(sec/2 + 1)), int(2*(1+N*sec/2)-2), -2):
        fk = np.delete(fk,i,0)
    for i in range(2*N, -2, -2):
        fk = np.delete(fk,i,0)
    return fk

## Calculate displacement
def cal_U(K, F):
    invK = np.linalg.inv(K)
    u = invK.dot(F)
    
    return u

## Make the dimension of displacement equal to the dimension of original F matrix
def Uremake(u, N, sec):
    z = np.zeros([1,1])
    U = np.insert(u,0,z,0)
    for i in range(2, 2*N +2, 2):
        U = np.insert(U,i,z,0)
    for i in range(int(2*(1+N*sec/2)), int(2*N*(sec/2 + 1) + 2),  2):
        U = np.insert(U,i,z,0)
    U = np.insert(U, int(2*N*(sec/2 + 1) + 1), z, 0)
    
    return U

## Determine coordinate of each node after deformation
def Deform(p, u):
    p = np.transpose(p)
    for i in range(len(p)):
        p[i][0] += u[2*i][0]

        if p[i][1] + u[2*i + 1][0] < -1:
            p[i][1] = -1
        else:
            p[i][1] += u[2*i + 1][0]

    return np.transpose(p)

TK = Ktot_pstrain(M, p , E, v )         # Total stiffness matrix
RK = Kreduce(TK, N, sec)                # Reduced K matrix
RF = Freduce(F, N, sec)                 # Reduced F matrix
u = cal_U(RK, RF)                       # Reduced displacement matrix
uu = Uremake(u, N, sec)                 # Displacement matrix
pp = Deform(p,uu)                       # Coordinate of all nodes after deformation

## Class for calculate stress under plane strain condition
class CalStress:
    def __init__(self, M, p, u, E, v):
        self.E = E
        self.v = v
        self.p = p
        self.M = M
        self.u = u

    def pstrain(self):
        
        sig = np.zeros([len(self.M), len(self.M[0])])
        
        # Making stiffness matrix
        for i in range(len(self.M)):
            D = K_pstrain(self.M[i], self.p, self.E, self.v).D()
            B = K_pstrain(self.M[i], self.p, self.E, self.v).B()

            ui = np.zeros([2*len(self.M[i]), 1])
            for j in range(len(self.M[i])):
                ui[2*j] = self.u[2*self.M[i][j]]
                ui[2*j +1] = self.u[2*self.M[i][j]+1]

            sig[i] = np.transpose(D.dot(B.dot(ui)))
        
        return sig
    
    def pstress(self):
        
        sig = np.zeros([len(self.M), len(self.M[0])])
        
        # Making stiffness matrix
        for i in range(len(self.M)):
            D = K_pstress(self.M[i], self.p, self.E, self.v).D()
            B = K_pstress(self.M[i], self.p, self.E, self.v).B()

            ui = np.zeros([2*len(self.M[i]), 1])
            for j in range(len(self.M[i])):
                ui[2*j] = self.u[2*self.M[i][j]]
                ui[2*j +1] = self.u[2*self.M[i][j]+1]

            sig[i] = np.transpose(D.dot(B.dot(ui)))
        
        return sig

    def Maxpstrain(self):
        Maxsig = np.zeros([2,3])
        for i in range(3):
            sigma = np.transpose(self.pstrain())[i]
            Maxsig[0][i] = max(sigma, key=abs)
            Maxsig[1][i] = np.where(sigma == max(sigma, key=abs))[0][0]

        return Maxsig

    def Maxpstress(self):
        Maxsig = np.zeros([2,3])
        for i in range(3):
            sigma = np.transpose(self.pstress())[i]
            Maxsig[0][i] = max(sigma, key=abs)
            Maxsig[1][i] = np.where(sigma == max(sigma, key=abs))[0][0]

        return Maxsig

## Stress matrix for all elements
sigma = CalStress(M, p, uu, E, v).pstrain()

sigmaXX = np.transpose(sigma)[0]
sigmaYY = np.transpose(sigma)[1]
sigmaXY = np.transpose(sigma)[2]

msh.Mesh(pp, M, sigmaXX, sigmaYY, sigmaXY)

