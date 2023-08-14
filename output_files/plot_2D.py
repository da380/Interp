import numpy as np
import matplotlib.pyplot as plt
import sys


#######################################################################

file = sys.argv[1]
xcs = sys.argv[2]
ycs = sys.argv[3]
file2 = sys.argv[4]
# xmv = sys.argv[4]
# xlab = sys.argv[5]
# ylab = sys.argv[6]
data = np.loadtxt(file, delimiter=";")
data2 = np.loadtxt(file2, delimiter=";")
xc = int(xcs)
yc = int(ycs)
# xmvc = float(xmv)
x = data[:, xc - 1]
y = data[:, yc - 1]
x2 = data2[:, 0]
y2 = data2[:, 1]
plt.plot(x, y, "k", linewidth=1.0)
plt.plot(x2, y2, "ro")
plt.show()
