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
from AC_Modules import *
warnings.filterwarnings(action = "ignore")
pd.set_option('display.max_rows', 10000)

# ===================================== MAIN BODY ==================================== #

if __name__ == "__main__": 
    
    ###############################################################################
    #                               DATA PROCESSING                               #  
    ###############################################################################

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
    k1 = rate_constant(rxn, 19, TEMP)
    J = rate_constant(rxn, 28, TEMP)
    k2 = rate_constant(rxn, 1, TEMP)
    
    ###############################################################################
    #                               INPUT/ITERATIONS                              #  
    ###############################################################################
    h_q, h_f, h_b = 0.1, 0.00001, 0.1 #time increment
    h = [h_q, h_f, h_b]
    T = 1000 # Total length of time
    
    
    ### ITERATIONS
    res_QSSA = QSSA(info, h_q, T, k1, J, k2)
    res_QSSA['NO+NO2'] = res_QSSA['NO'] + res_QSSA['NO2']
    res_FORWARD = FORWARD(info, h_f, T, k1, J, k2)
    res_FORWARD['NO+NO2'] = res_FORWARD['NO'] + res_FORWARD['NO2']
    res_BACKWARD = BACKWARD(info, h_b, T, k1, J, k2)
    res_BACKWARD['NO+NO2'] = res_BACKWARD['NO'] + res_BACKWARD['NO2']
        
    var_list = info[info['A/I'] == 'A']['NAME'].values

    t_axis_f = np.arange(0, T+h_f, h_f)
    t_axis_b = np.arange(0, T+h_b, h_b)
    t_axis_q = np.arange(0, T+h_q, h_q)

    t_axis = [t_axis_q, t_axis_f, t_axis_b]
    
    data_list = [res_QSSA, res_FORWARD, res_BACKWARD]
    method_text = ['QSSA', 'FORWARD', 'BACKWARD']
    data_text = ['QSSA_Diff.', 'FORWARD_Diff.', 'BACKWARD_Diff.']

    N_initial = res_QSSA['NO'][0] + res_QSSA['NO2'][0]
    

    clist = ['C1', 'C3', 'C0', 'C8']
    vlist = ['$NO$', '$NO_2$', '$O$', '$O_3$']


    
    ###############################################################################
    #                             MASS CONSERVATION                               #  
    ###############################################################################
    
    fig_mass, ax_mass = plt.subplots(1, 3, figsize = (20, 6))
    
    #Y_mass_max = [max(df['NO+NO2'].values) for df in data_list]
    #Y_mass_min = [min(df['NO+NO2'].values) for df in data_list]
    

    for i, df in enumerate(data_list):
        
        N_i = [N_initial]*round((T+h[i])/h[i])
        ax_mass[i].set_title(data_text[i])
        #ax_mass[i].plot(t_axis[i], N_i, label = 'Initial', color = 'C8', linestyle = '--')
        ax_mass[i].plot(t_axis[i], df['NO+NO2'] - N_i, color ='C0', label = data_text[i])
        ax_mass[i].plot(t_axis[i], np.zeros(len(t_axis[i])), color = 'C8', linestyle = '--')
        ax_mass[i].set_ylabel('[$NO$]+[$NO_2$], ($molecules$ $cm^{-3}$)')
        ax_mass[i].set_xlabel('Time(s)')
        ax_mass[i].legend(loc = 'lower right')
        ax_mass[i].set_xlim(0, T)
        #ax_mass[i].set_ylim(Y_mass_min[i], Y_mass_max[i])
        ax_mass[i].ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset = False)

    fig_mass.suptitle('Mass Conservation', fontweight = 'semibold', fontsize = 16, y = 0.98, x = 0.51)
    fig_mass.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig('Mass Conservation.png', dpi = 300)
    
    
    ###############################################################################
    #                                SUMMARY PLOTS                                #  
    ###############################################################################
    
    comp_list = np.append(var_list, res_BACKWARD.columns.values[-1])
    comp_text = ['$NO$', '$NO_2$', '$O$', '$O_3$', '$NO+NO_2$']
    ls_list = ['-','--',':']

    fig_comp, ax_comp = plt.subplots(1, 5, figsize = (30, 5))
    fig_comp.suptitle(' vs. '.join(method_text), fontweight = 'semibold', fontsize = 14)

    for i, var in enumerate(comp_list):
        
        ax_comp[i].set_title(comp_text[i])

        for j, df in enumerate(data_list):
            N_i = [N_initial]*round((T+h[j])/h[j])
            ax_comp[i].plot(t_axis[j], df[var], label = method_text[j], linestyle = ls_list[j])

        ax_comp[i].set_xlim(0, T)
        ax_comp[i].set_xlabel('Time[s]')
        ax_comp[i].set_ylabel(f'[{comp_text[i]}], $molecules\, cm^{-3}$')
        ax_comp[i].legend(loc = 'lower right')

    ax_comp[-1].plot(t_axis[0], [N_initial]*len(t_axis[0]), label = 'Initial', color = 'C8', linestyle = '-.')
    ax_comp[-1].legend(loc = 'lower right')

    plt.savefig('Final Comparison.png', dpi = 300)
    
    
    
    ###############################################################################
    #                               PLOT BY METHOD                                #  
    ###############################################################################
    for j, df in enumerate(data_list):
        fig, ax = plt.subplots(1,4, figsize = (22, 5))
        fig.suptitle(method_text[j], fontweight = 'semibold', x = 0.51)
        
        for i, var in enumerate(var_list):
            ax[i].plot(t_axis[j], df[var], label = vlist[i], color = clist[i] )
            ax[i].set_title(vlist[i])
            ax[i].set_ylabel(f'[{vlist[i]}]'+', $molecules\, cm^{-3}$')
            ax[i].set_xlabel(f'Time(s)\ntime step = {h[j]}')
            ax[i].legend(loc = 'lower right')
            
        fig.tight_layout()
        plt.savefig(method_text[j], dpi = 300)
        
  

