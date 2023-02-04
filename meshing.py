import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

pi = np.pi
R = 1
N = 4
theta = pi/4
r = R/N
mps = 2*N-1
sec = int(2*pi/theta)

p = np.zeros([2, N*sec + 1])
mesh = np.zeros([(mps)*sec, 3])

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
for i in range (0,sec):
    for j in range (1, N+1):
        p[0][N*i + j], p[1][N*i+j] = rotate([0.0, 0.0], [0.0, (j) * R/N], i * theta)
    
# print(np.shape(p))
# print(np.around(np.transpose(p),2))
# print(np.shape(mesh))

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

print(mesh)
print(np.shape(mesh))


# x = np.linspace(0, 1, 100)
# y = np.sqrt(1.0**2 - (x)**2)
plt.scatter(p[0], p[1])
triangulation = tri.Triangulation(p[0], p[1], mesh)
plt.triplot(triangulation, '-k')

# # plt.plot(p)
# plt.plot(-x, y ,'b')
# plt.plot(-x, -y ,'b')
plt.gca().set_aspect('equal')
plt.show()