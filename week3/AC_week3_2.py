import numpy as np
import matplotlib.pyplot as plt











# Contour Graph 

plt.figure(figsize = (7, 3))

#Subplot1
plt.subplot(121)
n1 = 100
x1 = np.linspace(-2*np.pi,2*np.pi,n1)
y1 = np.linspace(-2*np.pi,2*np.pi,n1)
X1, Y1 = np.meshgrid(x1,y1)

Z1 = np.sin(X1) + np.sin(Y1)
 
C1 = plt.contour(X1,Y1, Z1, cmap=plt.cm.Reds,linewidth = 0.5)
plt.clabel(C1, inline = 1, fontsize = 6)
plt.title('$ Z = sinX + sinY $')
plt.xlabel('X')
plt.ylabel('Y')

#Subplot2
plt.subplot(122)
n2 = 100
x2 = np.linspace(-2,2,n2)
y2 = np.linspace(-2,2,n2)
X2, Y2 = np.meshgrid(x2,y2)

Z2 = X2*np.exp(-X2**2-Y2**2)
C2 = plt.contour(X2, Y2, Z2, cmap = plt.cm.Blues,linewidth = 0.5)
plt.clabel(C2, inline = 1, fontsize = 6)
plt.title('$Z = Xe^{-X^2-Y^2}$')
plt.xlabel('X')
plt.ylabel('Y')








plt.show()



