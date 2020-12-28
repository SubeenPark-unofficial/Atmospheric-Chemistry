# =============================== GENERAL INFORMATIONS ================================= #
# Created by: Subeen Park
# Date: May 27, 2020
# Explanation: Numerical Computing Practice for Photostationary Ozone Chemistry (Forward + Backward + QSSA)

# =============================== INVOLVED REACTIONS ================================== #
# NOTATIONS: oxygen : O2 | ozone: O3 | oxygen radical : O3P | nitrogen oxide : NO | nitrogen dioxide : NO2
# rxn1: NO + O3 -> NO2 + O2 | rxn NMBR 19
# rxn2: NO2 + hv -> NO + O | rxn NMBR 28
# rxn3: O + O2 + M -> O3 + M | rxn NMBR 1

# =================================== IMPORT MODULES ================================== #
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numba import jit 
import warnings
from tqdm.notebook import tqdm
warnings.filterwarnings(action = "ignore")
pd.set_option('display.max_rows', 10000)


# ============================= FUNCTION DECLARATION(SKIP!) ============================== #

################################################################################
#                                 TIME FORMATTING                              #  
#####################################################c##########################

def sec_to_DHS(time_in_sec):
    d = time_in_sec//(3600*24)
    h = (time_in_sec - d*3600*24)//3600
    m = (time_in_sec - d*3600*24-h*3600)//60
    s = time_in_sec - d*(3600*24) - h*3600 -m*60
    return d, h, m, s


################################################################################
#                               CHEMICAL CONSTATNS                             #  
#####################################################c##########################
TEMP = 298
A_v = 6.022*10**23
P = 10**5 # unit: Pa = N/m^2
V = 10**(-6) # unit: m^2
T = TEMP # K
R = 8.31 # J/(molK)
conv = A_v*P*V/(R*T) # number of molecules of air in cubic centimiter

def rate_constant(df, rxn_num, T):
    line = df[df['NMBR'] == str(rxn_num)]
    A, B, C = float(line['A'].item()), float(line['B'].item()), float(line['C'].item())
    k = A * (300/T)**B * np.exp(C/T)
    return k

def diurnal_J(time_arr, J):
    J_arr = []
    A = 2*np.pi/86400
    for t in time_arr:
        hr = sec_to_DHS(t)[1]
        if hr>=6 and hr<= 17:
            J_arr.append(-J*np.cos(A*t))
        else:
            J_arr.append(0)
            
    return J_arr

###############################################################################
#                            1 TIME STEP INCREMENT                            #  
###############################################################################
@jit
def forward_nextstep(conc_t, P, L, h): 
    # conc_t = [X]_(t) , conc = [X]_(t+h) , P: Production, L = Loss, h = timestep
    conc = conc_t + h*(P-L)
    return conc

@jit
def backward_nextstep(conc_t, P, IL, h): 
    # conc_t = [X]_(t-h) , conc = [X]_(t) , P: Production, L =IMPLICIT LOSS!! , h = timestep
    conc = (conc_t + h*P)/(1+h*IL)
    return conc

@jit
def exp_nextstep(conc_t, P, IL, h):
    # conc_t = [X]_(t-h) , conc = [X]_(t) , P: Production, L =IMPLICIT LOSS!! , h = timestep
    conc = conc_t*np.exp(-h*IL) + (P/IL)*(1-np.exp(-h*IL))
    return conc

@jit
def steady_state(conc_t, P, IL):
    # conc_t = [X]_(t-h) , conc = [X]_(t) , P: Production, L =IMPLICIT LOSS!! , h = timestep
    conc = P/IL
    return conc


###############################################################################
#                            METHOD DETERMINATION                             #  
###############################################################################
def select_method(h, IL):
    
    if h*IL < 0.01:
        method = "forward"
    elif h*IL >= 0.01 and h*IL <= 10:
        method = "exponential"
    elif h*IL > 10:
        method = "steady-state"
        
    return method
    
def nextstep(conc_t, P, L, IL, h): 
    
    method = select_method(h, IL)
    
    if method == "forward":
        conc = forward_nextstep(conc_t, P, L, h)
    elif method == "exponential":
        conc = exp_nextstep(conc_t, P, IL, h)
    elif method == "steady-state":
        conc = steady_state(conc_t, P, IL)
        
    return conc

###############################################################################
#                               FORWARD METHOD                                #  
###############################################################################
@jit
def FORWARD(df, h, T, k1, J, k2):
    # df: information , T: lengh of time
    nt = round(T/h)
    
    # Create dataframe for calculation / Initialise data
    res = np.zeros((nt+1, 6))
    initial = conv*df['INITIAL'].values.astype(float)
    res[0,:] = initial

    
    # Set concentration of inactive species
    res[:,4].fill(initial[4]) 
    res[:,5].fill(initial[5])

    print ("START FORWARD....")        
    # Calculate Concentrations
    for step in range(1, nt + 1):
        
        # Concentration at previous step(Concentration at t)
        NO, NO2, O3P, O3, O2, M = res[step-1, :]
        
        # [NO]: Concentration at t+h
        P_NO = J*NO2
        L_NO = k1*NO*O3
        res[step,0] = forward_nextstep(NO, P_NO, L_NO, h)

        # [NO2]: Concentration at t+h
        P_NO2 = k1*NO*O3
        L_NO2 = J*NO2
        res[step,1] = forward_nextstep(NO2, P_NO2, L_NO2, h)

        # [O3P]: Concentration at t+h
        P_O3P = J*NO2
        L_O3P = k2*O3P*O2*M
        res[step,2] = forward_nextstep(O3P, P_O3P, L_O3P, h)
        
        # [O3]: Concentration at t+h
        P_O3 = k2*O3P*O2*M
        L_O3 = k1*NO*O3
        res[step,3] = forward_nextstep(O3, P_O3, L_O3, h)
    
    res = pd.DataFrame(res, columns = df['NAME'].values)
    print ( "!!!! CALCULATION IS DONE!!!!")
    print ("FORWARD METHOD:", T,  "seconds with dt =", h)
    return res


###############################################################################
#                               BACKWARD METHOD                               #  
###############################################################################
@jit
def BACKWARD(df, h, T, k1 , J , k2):
    # df: information , T: lengh of time
    nt = round(T/h)
    
    # Create dataframe for calculation / Initialise data
    res = np.zeros((nt+1, 6))
    initial = conv*df['INITIAL'].values.astype(float)
    res[0,:] = initial

    
    # Set concentration of inactive species
    res[:,4].fill(initial[4]) 
    res[:,5].fill(initial[5])
    
    print ("START BACKWARD....") 
    # Calculate Concentrations
    for step in range(1, nt + 1):
        
        # Concentration at previous step(Concentration at t-h)
        NO, NO2, O3P, O3, O2, M = res[step-1, :]
        
        # [NO]: Concentration at t
        P_NO = J*NO2
        L_NO = k1*O3
        res[step,0] = backward_nextstep(NO, P_NO, L_NO, h)

        # [NO2]: Concentration at t
        P_NO2 = k1*NO*O3
        L_NO2 = J
        res[step,1] = backward_nextstep(NO2, P_NO2, L_NO2, h)

        # [O3P]: Concentration at t
        P_O3P = J*NO2
        L_O3P = k2*O2*M
        res[step,2] = backward_nextstep(O3P, P_O3P, L_O3P, h)
        
        # [O3]: Concentration at t
        P_O3 = k2*O3P*O2*M
        L_O3 = k1*NO
        res[step,3] = backward_nextstep(O3, P_O3, L_O3, h)
        
    res = pd.DataFrame(res, columns = df['NAME'].values)
    print ( "!!!! CALCULATION IS DONE!!!!")
    print ("BACKWARD METHOD:", T,  "seconds with dt =", h)

    return res

###############################################################################
#                                  QSSA METHOD                                #  
###############################################################################
@jit
def QSSA(df, h, T, k1, J, k2):
    
    # df: information , T: lengh of time
    nt = round(T/h)
    
    # Create dataframe for calculation / Initialise data
    res = np.zeros((nt+1, 6))
    initial = conv*df['INITIAL'].values.astype(float)
    res[0,:] = initial


    # Set concentration of inactive species
    res[:,4].fill(initial[4]) 
    res[:,5].fill(initial[5])

    print ("START QSSA....") 
    # Calculate Concentrations
    for step in range(1, nt + 1):
                      
        NO, NO2, O3P, O3, O2, M = res[step-1, :]
        
    
        P_NO, L_NO, IL_NO = J*NO2, k1*NO*O3, k1*O3
        P_NO2, L_NO2, IL_NO2 = k1*NO*O3, J*NO2, J
        P_O3P, L_O3P, IL_O3P = J*NO2, k2*O3P*O2*M, k2*O2*M
        P_O3, L_O3, IL_O3 = k2*O3P*O2*M, k1*NO*O3, k1*NO
        
       
        res[step,0] = nextstep(NO, P_NO, L_NO, IL_NO, h)
        res[step,1] = nextstep(NO2, P_NO2, L_NO2, IL_NO2, h)
        res[step,2] = nextstep(O3P, P_O3P, L_O3P, IL_O3P, h)
        res[step,3] = nextstep(O3, P_O3, L_O3, IL_O3, h)
        
        
        
  
    res = pd.DataFrame(res, columns = df['NAME'].values)
    print ( "!!!! CALCULATION IS DONE!!!!")
    print ("QSSA METHOD:", T,  "seconds with dt =", h)
    
    return res

