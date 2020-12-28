import numpy as np
import matplotlib.pyplot as plt

# Saturated Vapour Pressure 
Temp_C = np.arange(1, 100, 0.5)
Temp_K = Temp_C - 273.15

es_1 = 6.11*np.exp(19.83 - 5417/Temp_K)
es_2 = 6.11*np.exp(17.625*Temp_C/(Temp_C+243.04))

plt.plot(Temp_C, es_2, 'bo-', label = 'Equation1')
plt.plot(Temp_C, es_2,'r--',label = 'Equation2')
plt.title("SVP vs. Temp")
plt.xlabel("Temperature[$^\circ$C]")
plt.ylabel("Saturated Vapour Pressure[mb]")
plt.legend(loc = 2)
plt.show()


