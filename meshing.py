import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

W = 1
H = 1

Nx = 3
Ny = 3

dx = W / (Nx-1)
dy = H / (Ny-1)

p = np.zeros([2,Nx*Ny])
mesh = np.zeros([2*(Nx-1)*(Ny-1), 3])

index = 0
for i in range (0,Ny):
    y = i*dy
    for j in range (-Nx+1,1):
        x = j*dx
        p[0][index] = x
        p[1][index] = y
        index+=1


print(p)

index = 0
for i in range (1,Ny):
    for j in range (1,Nx):
        index1 = j + (i-1)* Nx
        index2 = index1 + 1 
        index3 = index2 + Nx
        index4 = index1 + Nx
        
        mesh[index] = [index1, index3, index4]
        index += 1
        mesh[index] = [index1, index2, index3]
        index += 1 

print(mesh-1)

x = np.linspace(0, 1, 100)
y = np.sqrt(1.0**2 - (x)**2)
plt.scatter(p[0], p[1])
triangulation = tri.Triangulation(p[0], p[1], mesh-1)
plt.triplot(triangulation, '-k')

# plt.plot(p)
plt.plot(-x, y ,'b')
plt.plot(-x, -y ,'b')
plt.gca().set_aspect('equal')
plt.show()