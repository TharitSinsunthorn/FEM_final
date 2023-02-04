import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri


R = 1  # Radius of circle
N = 4  # number of division in r direction
theta = np.pi/4  # number of division in theta direction
mps = 2*N-1  # number of mesh per anglular division 
sec = int(2*np.pi/theta)  # number of angular division

p = np.zeros([2, N*sec + 1])  # Define matrix for coordinate
mesh = np.zeros([(mps)*sec, 3])  # Define matrix for mesh

def Findcoor(p, R, N, theta):
    def rotate(o, p, theta):
        """
        Rotate a point counterclockwise by a given angle around a given origin.

        The angle should be given in radians.
        """
        ox, oy = o
        px, py = p

        npx = ox + np.cos(theta) * (px - ox) - np.sin(theta) * (py - oy)
        npy = oy + np.sin(theta) * (px - ox) + np.cos(theta) * (py - oy)
        return npx, npy
    
    #Finding coordinate of all nodes
    for i in range (0, int(2*np.pi/theta)):
        for j in range (1, N+1):
            p[0][N*i + j], p[1][N*i+j] = rotate([0.0, 0.0], [0.0, (j) * R/N], i * theta)


def FindIndex(mesh, N, theta):
    sec = int(2*np.pi/theta)
    mps = 2*N-1
    # Finding index
    # mesh around origin
    mesh[mps*(sec-1)] = [0, 1 + N*(sec-1), 1]
    for i in range(0,sec-1):
        mesh[mps*i][0] = 0.0
        mesh[mps*i][1] = 1 + N*i
        mesh[mps*i][2] = 1 + N*(i+1)
    # mesh type1
    for i in range (0, sec):
        for j in range(1,N):
            m = [N*i + 1 + (j-1), N*i + 2 + (j-1), N*i + N+1 + (j-1)]

            for e in range(len(m)):
                if m[e] > N*sec:
                    m[e] = (m[e] - N*sec)
            
            mesh[i*(mps) + j] = m
    # mesh type2
    for i in range (0, sec):
        for j in range(1,N):
            m = [N*i + 2 + (j-1), N*(i+1) + 2 + (j-1), N*(i+1) + 1 + (j-1)]

            for e in range(len(m)):
                if m[e] > N*sec:
                    m[e] = (m[e] - N*sec)

            mesh[i*(mps) + j + N-1] = m



Findcoor(p, R, N, theta)
FindIndex(mesh, N, theta)
# print(mesh)

def Mesh(p, mesh):
    x = np.linspace(0, 1, 100)
    y = np.sqrt(1.0**2 - (x)**2)
    plt.scatter(p[0], p[1])
    triangulation = tri.Triangulation(p[0], p[1], mesh)
    plt.triplot(triangulation, '-k')

    # plt.plot(p)
    plt.plot(x, y ,'b')
    plt.plot(x, -y ,'b')
    plt.plot(-x, y ,'b')
    plt.plot(-x, -y ,'b')
    plt.gca().set_aspect('equal')
    plt.show()

# Mesh(p, mesh)