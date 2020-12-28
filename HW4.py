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
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, FuncFormatter)
from numba import jit 
from AC_Modules import *
#warnings.filterwarnings(action = "ignore")
pd.set_option('display.max_rows', 10000)

def diurnal_variation(df, h, T, k1, J_max, k2):
    nt = round(T/h)

    time_arr = np.arange(0, T+h, h)
    J_arr = diurnal_J(time_arr, J_max)

    # Create dataframe for calculation / Initialise data
    res = np.zeros((nt+1, 6))
    initial = conv*df['INITIAL'].values.astype(float)
    res[0,:] = initial


    # Set concentration of inactive species
    res[:,4].fill(initial[4]) 
    res[:,5].fill(initial[5])

    # Calculate Concentrations
    for step in range(1, nt + 1):

        # Concentration at previous step(Concentration at t-h)
        NO, NO2, O3P, O3, O2, M = res[step-1, :]
        J = J_arr[step]

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
    
    return res

# ===================================== MAIN BODY ==================================== #

if __name__ == "__main__":
    
    ### READ FILES
    data = pd.read_table('photo_stationary.dat', comment = '#', header = None)

    ### SPLIT DATA
    begins = data[data[0]=='BEGIN'].index.values.astype(int)
    ends = data[data[0]=='END'].index.values.astype(int)

    ### NAME/MW/BKGAS(VMRAT) 
    info = data.iloc[begins[0]+1:ends[0]][0]\
                .str.split(' +', expand = True)\
                .rename(columns = {0:'A/I', 1:'NAME', 2:'MW', 3:'INITIAL'})\
                .reset_index(drop = True)

    #### CHEMICAL REACTIONS
    chem = data.iloc[begins[1]+1:ends[1]][0].str.split(' +', expand = True).drop(columns = 7)
    chem = chem[chem[0]=='A'].rename(columns = {0:'Active', 1:'NMBR', 2:'A', 3:'B', 4:'C', 5:'Q', 6:'Fc'})

    photo = data.iloc[begins[2]+1:ends[2]][0].str.split(' +', expand = True)
    photo= photo[photo[0]=='A'].rename(columns = {0:'Active', 1:'NMBR', 2:'A', 3:'B', 4:'C', 5:'Q', 6:'Fc'})

    rxn = pd.concat([chem, photo]).reset_index(drop = True)
    
    ### RATE CONSTANTS
    TEMP = 298
    k1 = rate_constant(rxn, 19, TEMP)
    J_max = rate_constant(rxn, 28, TEMP)
    k2 = rate_constant(rxn, 1, TEMP)
    
    ### INPUT
    h = 1 #time increment
    T = 86400*3 # Total length of time

    ### ITERATIONS
    res_BACKWARD = diurnal_variation(info, h, T, k1 , J_max, k2)
    
    var_list = info['NAME'][info['A/I']=='A'].values
    var_text = ['$NO$', '$NO_2$', '$O$', '$O_3$']
    clist = ['C1', 'C3', 'C0', 'C8']
    t = np.arange(0, T+h, h)


    @FuncFormatter
    def major_formatter(x, pos):
        if sec_to_DHS(x)[0]<3:
            return f"Day{sec_to_DHS(x)[0]+1:.0f}"

    @FuncFormatter
    def minor_formatter(x, pos):
        return f"{int(sec_to_DHS(x)[1]):0>2}"

    fig, ax = plt.subplots(4,1, figsize = (8, 16))

    for i, var in enumerate(var_list):

        ax[i].plot(t, res_BACKWARD[var], color = clist[i])

        ax[i].set_xlim(0, T)
        ax[i].xaxis.set_major_locator(MultipleLocator(60*60*24))
        ax[i].xaxis.set_major_formatter(major_formatter)

        ax[i].xaxis.set_minor_locator(MultipleLocator(60*60*6))
        ax[i].xaxis.set_minor_formatter(minor_formatter)
        ax[i].tick_params('x', which = 'major', labelrotation = 0)
        ax[i].tick_params('x', which = 'minor', labelrotation = 0)

        ax[i].set_xlabel('Time(hr)')
        ax[i].set_ylabel(f'Concentration'+' ($molecules\, cm^{-3}$)')
        ax[i].set_title(var_text[i])
        ax[i].grid()


    fig.tight_layout()
    plt.savefig('HW4.png')
    plt.show()



