import numpy as np
import matplotlib.pyplot as plt
import sys


#######################################################################

file = sys.argv[1]

file2 = sys.argv[2]
# xmv = sys.argv[4]
# xlab = sys.argv[5]
# ylab = sys.argv[6]
data = np.loadtxt(file, delimiter=";")
data2 = np.loadtxt(file2, delimiter=";")

# xmvc = float(xmv)
x = data[:, 0]
y = data[:, 1]
y2 = data[:,2]
y3 = data[:,3]
x2 = data2[:, 0]
y4 = data2[:, 1]
plt.plot(x, y, "k", linewidth=1.0)
plt.plot(x, y2, "g", linewidth=1.0)
plt.plot(x, y3, "b", linewidth=1.0)
plt.plot(x2, y4, "ro")
plt.legend(["Exact", "Cubic", "Akima", "Interpolation Points"])
plt.show()
