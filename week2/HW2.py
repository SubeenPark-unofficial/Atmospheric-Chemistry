"""
Created on Thu May  7 12:13:27 2020

@author: subeen
"""


import numpy as np
import matplotlib.pyplot as plt

def psi(x):
    if abs(x-5/6) <= 1/9:
        return (9**4)*((x-5/6)**2 - (1/9)**2)**2
    else:
        return 0

def approx(X, T, dx, dt, U = 0.2):
    
    # Grid size 
    nx = int(X/dx)
    nt = int(T/dt)
    
    # X/psi array
    x = np.array([dx*i for i in range(nx+1)])
    t = np.array([dt*i for i in range(nt+1)])
    psi0 = np.array([psi(val) for val in x])
    

    # Initialization
    psi_LF = []
    psi_FW = []
    psi_LF.append(psi0)
    psi_FW.append(psi0)

    
    for n in range(0, nt):
        psi_t = []
        for i in range(0, nx+1):
            if i == nx:
                psi_t.append(psi_FW[n][i] - U*(dt/dx)*(psi_FW[n][1] - psi_FW[n][i]))
            else:
                psi_t.append(psi_FW[n][i] - U*(dt/dx)*(psi_FW[n][i+1] - psi_FW[n][i]))
                
        psi_FW.append(psi_t)
                              
    psi_LF.append(psi_FW[1])
    
    for n in range(1, nt):
        psi_t = []
        for i in range(0, nx+1):
            if i == 0:
                psi_t.append(psi_LF[n-1][i] - U*(dt/dx)*(psi_LF[n][i+1] - psi_LF[n][nx-1])) 
            elif i == nx:
                psi_t.append(psi_LF[n-1][i] - U*(dt/dx)*(psi_LF[n][1] - psi_LF[n][i-1]))
            else:
                psi_t.append(psi_LF[n-1][i] - U*(dt/dx)*(psi_LF[n][i+1] - psi_LF[n][i-1]))
        psi_LF.append(psi_t)
        
    plt.figure(figsize =(7, 5))    
    plt.plot(x, psi0, x, psi_LF[nt], x, psi_FW[nt])
    plt.xlim(0,1.0)
    plt.ylim(-0.6,1.6)
    plt.legend(['Analytic','Leapfrog', 'Forward'])
    plt.title(f"When dx = 1/{(1/dx):.0f}, dt = 1/{(1/dt):.0f}")
    plt.ylabel('$\Psi(x,5)$')
    plt.show()
    

approx(1, 5, 1/36, 1/72)
approx(1, 5, 1/36, 1/8)
approx(1, 5, 1/36, 1/12)