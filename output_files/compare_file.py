import numpy as np
import matplotlib.pyplot as plt
import sys


#######################################################################

file = sys.argv[1]
mylen = len(sys.argv)
# xmv = sys.argv[4]
# xlab = sys.argv[5]
# ylab = sys.argv[6]
data = np.loadtxt(file, delimiter=";")
# print (len(data[0,:]))
mylen = len(data[0,:])

# xmvc = float(xmv)
x = data[:, 0]
y = data[:, 1]
plt.plot(x, y, "k", linewidth=1.0)

if mylen >=3:
    y2 = data[:,2]
    plt.plot(x, y2, "g", linewidth=1.0)

if mylen >= 4:
    y3 = data[:,3]
    plt.plot(x, y3, "b", linewidth=1.0)

if mylen >= 5:
    y4 = data[:,4]
    plt.plot(x,y4,"y", linewidth = 1.0)



plt.legend(["Exact", "Cubic", "Akima"], fontsize = 20)
plt.show()
