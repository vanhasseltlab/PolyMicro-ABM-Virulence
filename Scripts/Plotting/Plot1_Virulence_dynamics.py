#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
An agent based model of a polymicrobial community with virulence-mediated interspecies interactions
Author: Eske van Meegen and Cathi Herzberg, LACDR, Leiden University, The Netherlands

Parts of the code are based on Herzberg C. (2023) Vanhasseltlab/polymicro-abm: Polymicro-abm. v1.0.0. An agent
based modeling framework of an interactive, polymicrobial community during antimicrobial treatment. Zenodo. 
doi:10.5281/zenodo.10039816 which can be found on https://github.com/vanhasseltlab/PolyMicro-ABM 

File Description: This script creates plots of single virulence factors in mono and dual-species populations. Figures 3, 4 and 5 (Filenames Virulence_Dynamics_['Ac', 'Ac + B_'].png, Virulence_Dynamics_['As', 'As + B_'].png, Virulence_Dynamics_['Ar', 'Ar + B_'].png)
"""
#%%

import os
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
#%%
#working directory: Main folder of repository "PolyMicro-ABM-Virulence"
path = os.getcwd() 
path_data = path + "/001_Complete_Database_Virulence_Strengths"
# load data
file = open(path_data + "/database_filled.pkl", "rb")
database = pickle.load(file)

#%% plot

# define plot formatting settings
plt.style.use('default')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['figure.dpi']= 600
mpl.rcParams['lines.linewidth']  = 0.5
mpl.rcParams['legend.fontsize'] = 5
mpl.rcParams['legend.title_fontsize'] = 7
mpl.rcParams['axes.labelsize']=7
mpl.rcParams['axes.titlesize']=7
mpl.rcParams['font.size']=7
mpl.rcParams['xtick.labelsize']=5 
mpl.rcParams['ytick.labelsize']=5
mpl.rcParams['axes.linewidth'] = 0.5
mpl.rcParams['xtick.major.size'] = 0.5
mpl.rcParams['xtick.major.width'] = 0.5
mpl.rcParams['xtick.minor.size'] = 0.5
mpl.rcParams['xtick.minor.width'] = 0.3
mpl.rcParams['ytick.major.size'] = 1
mpl.rcParams['ytick.major.width'] = 0.5
mpl.rcParams['ytick.minor.size'] = 0.5
mpl.rcParams['ytick.minor.width'] = 0.3
mpl.rcParams['grid.linewidth'] = 0.3
mpl.rcParams['lines.markersize'] = 0
mpl.rcParams['lines.marker'] = 'o'

# MAKE FIGURE FOLDER 
fig_folder = "/Figures/VirulenceDynamics" 
fig_path = path + fig_folder
if not (os.path.exists(fig_path)): # if the folder doesn't exist yet
    os.mkdir(fig_path)  

#%% plot each scenario individually
# colors and linestyles
mpl.rcParams['figure.figsize']=[3.2,3]

c=['red','orange','bisque','yellow','green','lime','deepskyblue','blue','violet','purple','fuchsia','darkred','orangered','greenyellow','lightgreen','cyan','navy','blueviolet','deeppink']
cmap = mpl.cm.get_cmap('Blues')
rgba1 = cmap(0.8)
cmap = mpl.cm.get_cmap('Oranges')
rgba2 = cmap(0.8)
c_species=[rgba1,rgba2,'grey']
runtime = 1440

vir_types=[['Ac','Ac + B_'],['As','As + B_'],['Ar','Ar + B_']]
v_pair=vir_types[0]
for v_pair in vir_types:
    fig = plt.figure()
    axes=[]
    outer = gridspec.GridSpec(1, 2, wspace=0.2, hspace=0)
    for key,i in zip(v_pair,range(2)):
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[i], wspace=0, hspace=0.5)
        
        #uppersubplot
        j=0
        x_alls=[database[key][key2][0]['mean']['Time'] for key2 in database[key]]
        
        y_alls=[database[key][key2][0]['mean']['M-Number_microbes'] for key2 in database[key]]
        ys_alls=[database[key][key2][0]['mean']['S-Number_microbes'] for key2 in database[key]]
        strengths=[key2 for key2 in database[key]]
        ax = plt.Subplot(fig, inner[j])
        axes=axes+[ax]
        c_i=0
        for x,y,y_s in zip(x_alls,y_alls,ys_alls):
            x=[x*10 for x in list(range(0,len(y)))]
            ax.plot(x,y,color=c[c_i])
            ax.fill_between(x,y-y_s,y+y_s,alpha=0.05,color=c[c_i])
            c_i+=1
        ax.set_xlim([0,runtime])
        ax.grid(which='major')
        ax.set_ylim([0,6000])
        ax.set_ylabel('Number of microbes')
        ax.set_xlabel('Time [min]')
        fig.add_subplot(ax)
        
        #lowersubplot
        j=1
        x=[float(key2[0]) for key2 in database[key]]
        
        ax = plt.Subplot(fig, inner[j])
        axes=axes+[ax]
        i=0
        for k in ['A','B']:
            if k in key:
                y=[database[key][key2][0]['summary_statistics']['mean']['m-AUC-#'+k] for key2 in database[key]]
                ax.plot(x,y,color=c_species[i],markersize=1)
                i+=1
        ax.set_ylim([0,7*10**6])
        ax.set_xlim([0,max([float(x[0]) for x in strengths])])
        if key=='As':
            ax.set_xlim([0,0.06])
        ax.set_ylabel('AUC')
        ax.set_xlabel('Virulence strength')
        ax.grid(which='major')
        fig.add_subplot(ax)
    
        
    axes[0].set_title('Mono-species \n A: '+'$v_{'+v_pair[0][1]+'}$')
    axes[2].set_title('Dual-species \n A: '+'$v_{'+v_pair[0][1]+'}$'+ '  +  B: -')
    axes[2].set_ylabel('')
    axes[3].set_ylabel('')
    axes[2].set_yticklabels([])
    axes[3].set_yticklabels([])
    y=1.2
    labels=['a)','b)','c)','d)']
    for i,l in zip([0,2,1,3],labels):
        ax=axes[i]
        axes[i].text(-0.15, y, l, fontsize=8,transform=ax.transAxes, verticalalignment='top')
    #add legend
    x=1.7
    legend_elements = []
    colors=c[0:len(strengths)]
    for s in range(len(strengths)):
        legend_elements += [Line2D([0], [0], color=colors[::-1][s], label=strengths[::-1][s][0])]
    ax=axes[3]
    ax.legend(handles=legend_elements,title='Virulence \n strengths',ncol=1,loc='lower right', bbox_to_anchor=(x, 1))

    legend_elements = []
    species=['virulent (A)','avirulent (B)']
    #species=['A','B']

    for s in range(len(species)):
        legend_elements += [Line2D([0], [0], color=c_species[s], label=species[s])]
    ax=axes[2]
    ax.legend(handles=legend_elements,title='Species',ncol=1,loc='lower right', bbox_to_anchor=(1.77, -1.3))
    

    plt.subplots_adjust(bottom=0.11, left=0.125, right=0.78, top=0.89)
        
    saveloc = fig_path + "/Virulence_Dynamics_"+str(v_pair)
    fig.savefig(saveloc+".png", format="png",dpi=600)