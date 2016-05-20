import gsw
import numpy as np
import matplotlib.pyplot as plt

#define S, T, p, lon, lat from your data

SA = gsw.SA_from_SP(S, p, lon, lat)
CT = gsw.CT_from_t(SA, T, p)

def SA_CT_plot(SA, CT, p_ref=0, isopycs=5, title_string=''):
    # if less than two input vars, error: "You need to supply both
    # Absolute Salinity and Conservative Temperature"
    # if len(p_ref)>1: error: Multiple reference pressures
    min_SA_data = np.amin(SA)
    max_SA_data = np.amax(SA)
    min_CT_data = np.amin(CT)
    max_CT_data = np.amax(CT)
    
    SA_min = max(0.0, min_SA_data-0.1*(max_SA_data - min_SA_data))
    SA_max = max_SA_data + 0.1*(max_SA_data - min_SA_data)
    SA_axis = np.arange(
        SA_min, SA_max + (SA_max-SA_min)/400, (SA_max-SA_min)/200)
    
    CT_freezing = gsw.CT_freezing(SA_axis, p_ref, 0)
    CT_min = min_CT_data - 0.1*(max_CT_data - min_CT_data)
    CT_max = max_CT_data + 0.1*(max_CT_data - min_CT_data)
    if CT_min > np.min(CT_freezing):
        CT_min = min_CT_data - 0.1*(max_CT_data - min(CT_freezing))
    CT_axis = np.arange(
        CT_min, CT_max + (CT_max-CT_min)/400, (CT_max-CT_min)/200)
    
    SA_gridded, CT_gridded = np.meshgrid(SA_axis, CT_axis)
    
    isopycs_gridded = gsw.rho(SA_gridded, CT_gridded,p_ref) - 1000.0
    
    c1 = plt.contour(
        SA_gridded, CT_gridded, isopycs_gridded, isopycs, colors='k')
    plt.clabel(c1, inline=1, fontsize=10)
    c2 = plt.plot(SA, CT, '.-', linewidth=2, markersize=10)
    
    # axis square?
    plt.axis((SA_min, SA_max, CT_min, CT_max))
    plt.xlabel('Absolute Salinity $\it{S}_A$ (g kg$^{-1}$)')
    plt.ylabel('Conservative Temperature, $\Theta$ ($^\circ$C)')
    if len(title_string) > 0:
        plt.title(title_string)
    else:
        plt.title('$\it{S}_A$ - $\Theta$ diagram:' + ' p$_{ref}$ = '
                  + str(p_ref) + ' dbar')
    plt.plot(SA_axis, CT_freezing, '--')
    #plt.gca().annotate('p$_{ref}$ = '+str(p_ref)+' dbar', xy=(2, 15), 
    #xytext=(2,20))
    
SA_CT_plot(SA, CT)
