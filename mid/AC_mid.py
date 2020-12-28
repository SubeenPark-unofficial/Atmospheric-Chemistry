# =============================== GENERAL INFORMATIONS ================================= #
# Created by: Subeen Park
# Date: May 20, 2020
# Explanation: Numerical Computing Practice for Photostationary Ozone Chemistry

# =============================== INVOLVED REACTIONS ================================== #
# NOTATIONS: oxygen : O2 | ozone: O3 | oxygen radical : O3P | nitrogen oxide : NO | nitrogen dioxide : NO2
# rxn1: NO + O3 -> NO2 + O2 | rxn NMBR 19
# rxn2: NO2 + hv -> NO + O | rxn NMBR 28
# rxn3: O + O2 + M -> O3 + M | rxn NMBR 1

# =================================== IMPORT MODULES ================================== #
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
pd.set_option('display.max_rows', 10000)

# ===================================== READ FILES ==================================== #
data = pd.read_table('photo_stationary.dat', comment = '#', header = None)

# =================================== SPLIT DATA ====================================== #
begins = data[data[0]=='BEGIN'].index.values.astype(int)
ends = data[data[0]=='END'].index.values.astype(int)

# NAME/MW/BKGAS(VMRAT) 

info = data.iloc[begins[0]+1:ends[0]][0]\
            .str.split(' +', expand = True)\
            .rename(columns = {0:'A/I', 1:'NAME', 2:'MW', 3:'INITIAL'})\
            .reset_index(drop = True)

# CHEMICAL REACTIONS
chem = data.iloc[begins[1]+1:ends[1]][0].str.split(' +', expand = True).drop(columns = 7)
chem = chem[chem[0]=='A'].rename(columns = {0:'Active', 1:'NMBR', 2:'A', 3:'B', 4:'C', 5:'Q', 6:'Fc'})

photo = data.iloc[begins[2]+1:ends[2]][0].str.split(' +', expand = True)
photo= photo[photo[0]=='A'].rename(columns = {0:'Active', 1:'NMBR', 2:'A', 3:'B', 4:'C', 5:'Q', 6:'Fc'})

rxn = pd.concat([chem, photo]).reset_index(drop = True)

# =================================== RATE CONSTANT ==================================== #
# Rate constants have form K = A * (300/T)**B * EXP(C/T)
TEMP = 298 #unit: K

def rate_constant(df, rxn_num, T = TEMP):
    line = df[df['NMBR'] == str(rxn_num)]
    A, B, C = float(line['A'].item()), float(line['B'].item()), float(line['C'].item())
    k = A * (300/T)**B * np.exp(C/T)
    return k

k1 = rate_constant(rxn, 19)
J = rate_constant(rxn, 28)
k2 = rate_constant(rxn, 1)

# ============================= UNIT CONVERSION: C_x to n_x =========================== #
A_v = 6.022*10**23
P = 10**5 # unit: Pa = N/m^2
V = 10**(-6) # unit: m^2
T = TEMP # K
R = 8.31 # J/(molK)

conv = A_v*P*V/(R*T) # number of molecules of air in cubic centimiter

    
# ================================= FORWARD METHOD ==================================== #
def forward_nextstep(conc_t, P, L, h): 
    # conc_t = [X]_(t) , conc = [X]_(t+h) , P: Production, L = Loss, h = timestep
    conc = conc_t + h*(P-L)
    return conc


def FORWARD(df, h, T):
    # df: information , T: lengh of time
    nt = int(T/h)
    
    # Create dataframe for calculation / Initialise data
    res = np.zeros((nt+1, 6))
    initial = conv*df['INITIAL'].values.astype(float)
    res[0,:] = initial

    
    # Set concentration of inactive species
    res[:,4].fill(initial[4]) 
    res[:,5].fill(initial[5])

            
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
   
    return res


# =============================== BACKWARD METHOD ================================== #
def backward_nextstep(conc_t, P, L, h): 
    # conc_t = [X]_(t-h) , conc = [X]_(t) , P: Production, L =IMPLICIT LOSS!! , h = timestep
    conc = (conc_t + h*P)/(1+h*L)
    return conc

def BACKWARD(df, h, T):
    # df: information , T: lengh of time
    nt = int(T/h)
    
    # Create dataframe for calculation / Initialise data
    res = pd.DataFrame(index = range(nt+1), columns = df['NAME'].values)
    res.loc[0] = conv*df['INITIAL'].values.astype(float)
    
    # Set concentration of inactive species
    for species in res.columns.values:
        if df[df['NAME'] == species]['A/I'].values[0] == 'I':
            res[species] = res[species].loc[0]
                    
    # Calculate Concentrations
    for step in range(1, nt + 1):
        
        # Concentration at previous step(Concentration at t-h)
        NO, NO2, O3P, O3, O2, M = res.loc[step-1]
        
        # [NO]: Concentration at t
        P_NO = J*NO2
        L_NO = k1*O3
        res['NO'].loc[step] = backward_nextstep(NO, P_NO, L_NO, h)

        # [NO2]: Concentration at t
        P_NO2 = k1*NO*O3
        L_NO2 = J
        res['NO2'].loc[step] = backward_nextstep(NO2, P_NO2, L_NO2, h)

        # [O3P]: Concentration at t
        P_O3P = J*NO2
        L_O3P = k2*O2*M
        res['O3P'].loc[step] = backward_nextstep(O3P, P_O3P, L_O3P, h)
        
        # [O3]: Concentration at t
        P_O3 = k2*O3P*O2*M
        L_O3 = k1*NO

        res['O3'].loc[step] = backward_nextstep(O3, P_O3, L_O3, h)
   
    return res

# ================================ ITERATION =================================== #
res_f = FORWARD(info, 10**(-5), 700)
res_b = BACKWARD(info, 0.01, 700)

# ================================ PLOT FIGURE ================================== #
t_f = np.arange(0,700 + 10**(-5),10**(-5))
O3_f = res_f['O3']

t_b = np.arange(0,700.01,0.01)
O3_b = res_b['O3']

O3_ss = 4.84*10**9*np.ones(len(t_b))

plt.figure(figsize = (8, 5))
plt.plot(t_b, O3_ss, 'b--', label = 'Theoretical')
plt.plot(t_b, O3_b, 'g', label = 'Backward')
plt.plot(t_f, O3_f, 'r:', label = 'Forward')
plt.title('Concentration of Ozone' , fontsize = 15)
plt.xlabel('time(s)\n timestep : h_forward = 0.00001s, h_backward = 0.01s')
plt.xlim(0,700)
plt.ylabel('$O_3$[$molecules \,cm^{-3}$]')
plt.legend(loc = 'lower right')

plt.savefig('MID_O3_highres.png', dpi = 300)
plt.show()
