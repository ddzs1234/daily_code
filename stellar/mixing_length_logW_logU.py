import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np

from sympy.solvers import solve
from sympy import sqrt,log
from sympy import Symbol

q=Symbol('q')
# make these smaller to increase the resolution
dx, dy = 1, 1

# generate 2 2d grids for the x & y bounds
y, x = np.mgrid[slice(-6, 6 + dy, dy),
                slice(-6, 6 + dx, dx)]

z=np.zeros(x.shape)

for i in range(0,13):
    for j in range(0,13):
#         print(i,j)
        x_ij=np.float_power(10,x[i,j])
        y_ij=np.float_power(10,y[i,j])
        m=solve((sqrt(q+y_ij**2)-y_ij)**3+8*y_ij*(q-x_ij)/9,q)[0]
        z[i,j]=log(m,10)
# print(x.shape,y.shape,z)
# # x and y are bounds, so z should be the value *inside* those bounds.
# # Therefore, remove the last value from the z array.
# z = z[:-1, :-1]
levels = MaxNLocator(nbins=50).tick_values(z.min(), z.max())


# # pick the desired colormap, sensible levels, and define a normalization
# # instance which takes data values and translates those into levels.
cmap = plt.get_cmap('plasma')


fig, (ax1) = plt.subplots(nrows=1)

# im = ax0.pcolormesh(x, y, z, cmap=cmap, norm=norm)
# fig.colorbar(im, ax=ax0)
# ax0.set_title('pcolormesh with levels')


# # contours are *point* based plots, so convert our bound into point
# # centers
cf = ax1.contourf(y,
                  x, z, levels=levels,
                  cmap=cmap)
fig.colorbar(cf, ax=ax1,label='logx')
ax1.set_xlabel('logU')
ax1.set_ylabel('logW')



# # adjust spacing between subplots so `ax1` title and `ax0` tick labels
# don't overlap
fig.tight_layout()

plt.show()
