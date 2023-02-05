import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.tri as tri


R = 1                            ## Radius of circle
N = 4                            ## number of division in r direction
theta = np.pi/4                  ## number of division in theta direction
mps = 2*N-1                      ## number of mesh per anglular division 
sec = int(2*np.pi/theta)         ## number of angular division

p = np.zeros([2, N*sec + 1])     ## Define matrix for coordinate
mesh = np.zeros([(mps)*sec, 3])  ## Define matrix for mesh


## Find coordinate of each node
def Findcoor(p, R, N, theta):
    def rotate(o, p, theta):
        ox, oy = o
        px, py = p

        npx = ox + np.cos(theta) * (px - ox) - np.sin(theta) * (py - oy)
        npy = oy + np.sin(theta) * (px - ox) + np.cos(theta) * (py - oy)
        return npx, npy
    
    for i in range (0, int(2*np.pi/theta)):
        for j in range (1, N+1):
            p[0][N*i + j], p[1][N*i+j] = rotate([0.0, 0.0], [0.0, (j) * R/N], i * theta)

## Define index for all nodes to each mesh
def FindIndex(mesh, N, theta):
    sec = int(2*np.pi/theta)
    mps = 2*N-1

    # Finding index
    # mesh around origin
    mesh[mps*(sec-1)] = [0, 1 + N*(sec-1), 1]
    for i in range(0,sec-1):
        mesh[mps*i][0] = 0
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

    mesh.astype(int)


Findcoor(p, R, N, theta)
FindIndex(mesh, N, theta)

## Simulation path
def Mesh(p, mesh, xx, yy, xy):
    
    ## Figure 1 shows contour of sigma XX
    Fxx = plt.figure(1)
    plt.title("Element: " + str(len(mesh)) + "         Node: " + str(len(p[0])))

    # plot circle
    x = np.linspace(-1, 1, 100)
    y = np.sqrt(1.0**2 - (x)**2)
    plt.plot(x, y ,'k')
    plt.plot(x, -y ,'k')

    # plot floor
    plt.axhline(y=-1, linewidth = 3, color = 'k', linestyle = '-')

    # Triangulation for meshing
    triangulation = tri.Triangulation(p[0], p[1], mesh)
    plt.triplot(triangulation, 'k-', marker = 'o', markerfacecolor = 'k')
    # Make a color contour for element's stress
    plt.tripcolor(p[0], p[1], mesh, alpha = 1.0, cmap = 'bwr', facecolors=xx, edgecolors='k')
    plt.colorbar(label = r'$\sigma_{xx}$')
    plt.gca().set_aspect('equal')

    ## Figure 2 shows contour of sigma YY
    Fyy = plt.figure(2)
    plt.title("Element: " + str(len(mesh)) + "         Node: " + str(len(p[0])))

    # plot circle
    x = np.linspace(-1, 1, 100)
    y = np.sqrt(1.0**2 - (x)**2)
    plt.plot(x, y ,'k')
    plt.plot(x, -y ,'k')

    # plot floor
    plt.axhline(y=-1, linewidth = 4, color = 'k', linestyle = '-')

    # Triangulation for meshing
    triangulation = tri.Triangulation(p[0], p[1], mesh)
    plt.triplot(triangulation, 'k-', marker = 'o', markerfacecolor = 'k')
    # Make a color contour for element's stress
    plt.tripcolor(p[0], p[1], mesh, alpha = 1.0, cmap = 'bwr', facecolors=yy, edgecolors='k')
    plt.colorbar(label = r'$\sigma_{yy}$')
    plt.gca().set_aspect('equal')

    ## Figure 3 shows contour of sigma XY
    Fxy = plt.figure(3)
    plt.title("Element: " + str(len(mesh)) + "         Node: " + str(len(p[0])))

    # plot circle
    x = np.linspace(-1, 1, 100)
    y = np.sqrt(1.0**2 - (x)**2)
    plt.plot(x, y ,'k')
    plt.plot(x, -y ,'k')

    # plot floor
    plt.axhline(y=-1, linewidth = 4, color = 'k', linestyle = '-')

    # Triangulation for meshing
    triangulation = tri.Triangulation(p[0], p[1], mesh)
    plt.triplot(triangulation, 'k-', marker = 'o', markerfacecolor = 'k')
    # Make a color contour for element's stress
    plt.tripcolor(p[0], p[1], mesh, alpha = 1.0, cmap = 'bwr', facecolors=xy, edgecolors='k')
    plt.colorbar(label = r'$\sigma_{xy}$')
    plt.gca().set_aspect('equal')    
   
    
    plt.show()


# Mesh(p, mesh)