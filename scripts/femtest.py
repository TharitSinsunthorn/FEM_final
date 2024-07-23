

# # A = np.array([[ 1166.67,     0.,   -1066.67,   133.33,     0.,    -333.33,  -100.,     200.  ],
# #  [    0.,     666.67,   200.,    -400. ,   -333.33,     0. ,    133.33,  -266.67],
# #  [-1066.67,   200. ,   1166.67,  -333.33 , -100.  ,   133.33  ,   0. ,      0.  ],
# #  [  133.33 , -400.  ,  -333.33 ,  666.67 ,  200. ,   -266.67  ,   0. ,      0.  ],
# #  [    0. ,   -333.33,  -100. ,    200. ,   1166.67 ,    0. ,  -1066.67,   133.33],
# #  [ -333.33 ,    0. ,    133.33 , -266.67 ,    0. ,    666.67 ,  200. ,   -400.  ],
# #  [ -100. ,    133.33 ,    0. ,      0. ,  -1066.67,   200. ,   1166.67 , -333.33],
# #  [  200. ,   -266.67  ,   0. ,      0.  ,   133.33 , -400. ,   -333.33 ,  666.67]])


import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import matplotlib.image as mpimg

# node_coordinate = {1: [0.0, 1.0], 2: [0.0, 0.0], 3: [4.018905, 0.87781],
#                    4: [3.978008, -0.1229], 5: [1.983549, -0.038322],
#                    6: [2.013683, 0.958586], 7: [3.018193, 0.922264],
#                    8: [2.979695, -0.079299], 9: [1.0070439, 0.989987],
#                    10: [0.9909098, -0.014787999999999999]}
# element_stress = {1: 0.2572e+01, 2: 0.8214e+00, 3: 0.5689e+01,
#                   4: -0.8214e+00, 5: -0.2572e+01, 6: -0.4292e+01,
#                   7: 0.4292e+01, 8: -0.5689e+01}

# n = len(element_stress.keys())
# x = np.empty(n)
# y = np.empty(n)
# d = np.empty(n)

# for i in element_stress.keys():
#     x[i-1] = node_coordinate[i][0]
#     y[i-1] = node_coordinate[i][1]
#     d[i-1] = element_stress[i]

# mask = np.logical_or(x < 1.e20, y < 1.e20)
# x = np.compress(mask, x)
# y = np.compress(mask, y)
# triang = tri.Triangulation(x, y)
# cmap = mpl.cm.jet
# fig = plt.figure(figsize=(80, 40))
# ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])

# levels=np.linspace(d.min(), d.max(), num=101)
# tri = ax1.tricontourf(triang, d, cmap=cmap, levels=levels)
# fig.colorbar(tri)
# plt.show()

h = 300
w = 1000
npts = 500
pts = np.zeros((npts,2))
pts[:,0] = np.random.randint(0,w,npts)
pts[:,1] = np.random.randint(0,h,npts)
tri = Delaunay(pts)
plt.xlim(0, w)
plt.ylim(0, h)
centers = np.sum(pts[tri.simplices], axis=1, dtype='int')/3.0
colors = np.array([ (x-w/2.)**2 + (y-h/2.)**2 for x,y in centers])
plt.tripcolor(pts[:,0], pts[:,1], tri.simplices.copy(), facecolors=colors, edgecolors='k')
plt.gca().set_aspect('equal')
# plt.show()
print(len(pts[:,0]))
print(len(tri.simplices.copy()))
print(len(colors))