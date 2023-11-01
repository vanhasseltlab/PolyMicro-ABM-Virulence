# -*- coding: utf-8 -*-
"""
An agent based model of a polymicrobial community with virulence-mediated interspecies interactions
Author: Eske van Meegen and Cathi Herzberg, LACDR, Leiden University, The Netherlands

Parts of the code are based on Herzberg C. (2023) Vanhasseltlab/polymicro-abm: Polymicro-abm. v1.0.0. An agent
based modeling framework of an interactive, polymicrobial community during antimicrobial treatment. Zenodo. 
doi:10.5281/zenodo.10039816 which can be found on https://github.com/vanhasseltlab/PolyMicro-ABM 

File Description: This script creates AUC heatmaps of dual-species populations. Figures 6 and 7
(Filenames AUC_v3s_[['Ac + B_', 'As + B_', 'Acs + B_'], ['Ac + Bs', 'Acs + Bc', 'Acs + Bs']].png, 
 AUC_v3s_[['As + B_', 'Asr + B_', 'As + Br'], ['Acs + B_', 'Acsr + B_', 'Acs + Br'], ['Asr + Br', 'Asr + Bs']].png)
"""
#%%
import os
import pickle
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#%%

# mpl.rcParams # check this dictionary for all the tweakable options
plt.style.use('default')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['figure.dpi']= 600
mpl.rcParams['figure.figsize']=[6.52,5.5]
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

cm = ['Blues', 'Oranges']

#%% OPEN COMPLETE DATABASE

#working directory: Main folder of repository "PolyMicro-ABM-Virulence"
path = os.getcwd() 
path_data = path + "/00_Complete_Database"
# load data
file = open(path_data + "/database_filled.pkl", "rb")
database = pickle.load(file)

fig_folder = "/Figures/AUCs_together" 
fig_path = path + fig_folder
if not (os.path.exists(fig_path)): # if the folder doesn't exist yet
    os.mkdir(fig_path)  
#%% PLOT: multiple heatmaps as subplots into one figure
column_names = ['m-AUC-#A', 'm-AUC-#B']
abconc = 0
species=['A','B']
mpl.rcParams['figure.figsize']=[6.52,3.3]

# get max AUC per species
all_aucs = [[],[]]
for virtype in database:
    for virstrength in database[virtype]:
        for s_i, column in enumerate(column_names):
            auc_i = database[virtype][virstrength][abconc]['summary_statistics']['mean'][column].loc[0] 
            all_aucs[s_i].append(auc_i)

vmax_AUCs = [[],[]]
for s_i, column in enumerate(column_names):
    vmax_AUCs[s_i] = max(all_aucs[s_i])

virtypes1 = [['Ac + B_', 'As + B_', 'Acs + B_'], ['Ac + Bs', 'Acs + Bc', 'Acs + Bs']]

virtypes2 = [['As + B_', 'Asr + B_', 'As + Br'], ['Acs + B_', 'Acsr + B_', 'Acs + Br'], ['Asr + Br', 'Asr + Bs']]
vs=[virtypes1,virtypes2]
vs=[virtypes1]    
for v in vs:
    fig = plt.figure()
    axes=[]
    
    outer = gridspec.GridSpec(len(v), 1, wspace=0, hspace=0.3) # create one subplot per row
    cbar_ax1 = fig.add_axes([.25, .12, .5, .06]) #add axis for scale bar
    cbar_ax2 = fig.add_axes([.25, .05, .5, .06]) #add axis for scale bar
    for v_row,i in zip(v,range(len(v))):
        inner = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[i], wspace=0.6, hspace=0) #create one subsubplot per item in row
        for virtype,k in zip(v_row,range(len(v_row))):
            inner_inner = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=inner[k], wspace=0.05, hspace=0) #create one subsubplot for each species
            for  j in range(2):
                s_i=j
                ax = plt.Subplot(fig, inner_inner[j])
                axes=axes+[ax]
                #collect correct auc data
                aucs = [[],[],[]]              
                for i, virstrength in enumerate(database[virtype]):
                    if len(database[virtype]) == 3: # B is avirulent: row,column = 1,3
                        x = 0
                    elif len(database[virtype]) == 9: # B is virulent: row,column = 3,3
                        x = i%3
                    else: 
                        raise Exception('length virtype not suitable for 1,3 or 3,3 colormap') 
                    auc_i = database[virtype][virstrength][abconc]['summary_statistics']['mean'][column_names[s_i]].loc[0]
                    aucs[x].append(auc_i)  
                    
                aucs = [ele for ele in aucs if ele != []]    # if row,column = 1,3 remove empty lists from aucs
                #flip order of rows
                aucs=[x for x in aucs[::-1]]
                #plot in heatmaps
                im = ax.imshow(np.array(aucs), vmin=0, vmax=vmax_AUCs[j], cmap=cm[j])
                if j==0:
                    cbar = plt.colorbar(im,cax=cbar_ax1, orientation='horizontal')
                    cbar_ax1.xaxis.set_ticks_position('top')
                if j==1:
                    cbar = plt.colorbar(im,cax=cbar_ax2, orientation='horizontal')
                    cbar_ax2.xaxis.set_ticks_position('bottom')
                vir_str = ' + '.join(['$v_{'+k+'}$' for k in virtype.split(' + ')[j][1::] if k!='_']) # creates v_ title                
                
                vir_x = ' + '.join(['$v_{'+k+'}$' for k in virtype.split(' + ')[0][1::] if k!='_']) # creates v_ title
                vir_y =  ' + '.join(['$v_{'+k+'}$' for k in virtype.split(' + ')[1][1::] if k!='_']) # creates v_ title

                title_string=species[j]+': '+' + '.join(['$v_{'+k+'}$' for k in virtype.split(' + ')[j][1::] if k!='_'])
                if title_string[-1]==' ':
                    title_string+='-'
                ax.set_title(title_string+'\n')
                if j==0:    
                    ax.set_xlabel('                     strength of A')
                    if len(vir_y)>0:    
                        ax.set_ylabel('strength of B')
                    else:
                        ax.set_yticks([])
                        # ax.set_ylabel('avirulent \n species B')
                fig.add_subplot(ax)
                
                plt.subplots_adjust(bottom=0.3, left=0.05, right=0.99, top=0.95)
                
                
    labels=['a)','b)','c)']
    for i,l in zip([0,2,4],labels):
        ax=axes[i]
        axes[i].text(-0.3, 2.8, l, fontsize=8,transform=ax.transAxes, verticalalignment='top')
    labels=['d)','e)','f)'] 
    for i,l in zip([6,8,10],labels):
        ax=axes[i]
        axes[i].text(-0.3, 1.5, l, fontsize=8,transform=ax.transAxes, verticalalignment='top')
           
    # remove y ticks on species B            
    for k in range(1,len(axes),2):
        axes[k].set_yticks([])
    
    # replace ytick labels with L, M, H
    for k in range(len(axes)):
        if len(axes[k].get_ylabel())>0:
            if axes[k].get_ylabel()[0]=='s':
                axis_ticks_y = ['$H$', '$M$', '$L$']
                axes[k].set_yticks(np.arange(len(axis_ticks_y)), labels=axis_ticks_y)
            else:
                axes[k].set_yticklabels([])
                
    # x tick labels
    for k in range(len(axes)):
        axis_ticks_x = ['$L$', '$M$', '$H$']
        axes[k].set_xticks(np.arange(len(axis_ticks_x)), labels=axis_ticks_x)

    props = dict(boxstyle='round', facecolor='white', alpha=0.1,edgecolor='white')
    axes[3].text(-6.0,-0.85, 'Normalized AUC \n of each species', transform=ax.transAxes, verticalalignment='top', bbox=props)
    axes[3].text(0,-0.8, 'A\n \nB', transform=ax.transAxes, verticalalignment='top', bbox=props)
    
    saveloc = fig_path + "/AUC_v3s_"+str(v)
    fig.savefig(saveloc+".png", format="png",dpi=600)
#%% PLOT: multiple heatmaps as subplots into one figure
column_names = ['m-AUC-#A', 'm-AUC-#B']
abconc = 0
species=['A','B']
mpl.rcParams['figure.figsize']=[6.52,5.5]
# get max AUC per species
all_aucs = [[],[]]
for virtype in database:
    for virstrength in database[virtype]:
        for s_i, column in enumerate(column_names):
            auc_i = database[virtype][virstrength][abconc]['summary_statistics']['mean'][column].loc[0] 
            all_aucs[s_i].append(auc_i)

vmax_AUCs = [[],[]]
for s_i, column in enumerate(column_names):
    vmax_AUCs[s_i] = max(all_aucs[s_i])

virtypes1 = [['Ac + B_', 'As + B_', 'Acs + B_'], ['Ac + Bs', 'Acs + Bc', 'Acs + Bs']]

virtypes2 = [['As + B_', 'Asr + B_', 'As + Br'], ['Acs + B_', 'Acsr + B_', 'Acs + Br'], ['Asr + Br', 'Asr + Bs']]
vs=[virtypes1,virtypes2]
vs=[virtypes2]    
for v in vs:
    fig = plt.figure()
    axes=[]
    
    outer = gridspec.GridSpec(len(v), 1, wspace=0, hspace=0.7) # create one subplot per row
    cbar_ax1 = fig.add_axes([.25, .10, .5, .04]) #add axis for scale bar
    cbar_ax2 = fig.add_axes([.25, .05, .5, .04]) #add axis for scale bar
    for v_row,i in zip(v,range(len(v))):
        inner = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[i], wspace=0.6, hspace=0) #create one subsubplot per item in row
        for virtype,k in zip(v_row,range(len(v_row))):
            inner_inner = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=inner[k], wspace=0.05, hspace=0) #create one subsubplot for each species
            for  j in range(2):
                s_i=j
                ax = plt.Subplot(fig, inner_inner[j])
                axes=axes+[ax]
                #collect correct auc data
                aucs = [[],[],[]]              
                for i, virstrength in enumerate(database[virtype]):
                    if len(database[virtype]) == 3: # B is avirulent: row,column = 1,3
                        x = 0
                    elif len(database[virtype]) == 9: # B is virulent: row,column = 3,3
                        x = i%3
                    else: 
                        raise Exception('length virtype not suitable for 1,3 or 3,3 colormap') 
                    auc_i = database[virtype][virstrength][abconc]['summary_statistics']['mean'][column_names[s_i]].loc[0]
                    aucs[x].append(auc_i)  
                    
                aucs = [ele for ele in aucs if ele != []]    # if row,column = 1,3 remove empty lists from aucs  
                # print(aucs)
                aucs=[x for x in aucs[::-1]]#flip rows
                #plot in heatmaps
                im = ax.imshow(np.array(aucs), vmin=0, vmax=vmax_AUCs[j], cmap=cm[j])
                if j==0:
                    cbar = plt.colorbar(im,cax=cbar_ax1, orientation='horizontal')
                    cbar_ax1.xaxis.set_ticks_position('top')
                if j==1:
                    cbar = plt.colorbar(im,cax=cbar_ax2, orientation='horizontal')
                    cbar_ax2.xaxis.set_ticks_position('bottom')
                vir_str = ' + '.join(['$v_{'+k+'}$' for k in virtype.split(' + ')[j][1::] if k!='_']) # creates v_ title                
                # ax.set_title(species[j])
                vir_x = ' + '.join(['$v_{'+k+'}$' for k in virtype.split(' + ')[0][1::] if k!='_']) # creates v_ title
                vir_y =  ' + '.join(['$v_{'+k+'}$' for k in virtype.split(' + ')[1][1::] if k!='_']) # creates v_ title
                title_string=species[j]+': '+' + '.join(['$v_{'+k+'}$' for k in virtype.split(' + ')[j][1::] if k!='_'])
                if title_string[-1]==' ':
                    title_string+='-'
                ax.set_title(title_string+'\n')
                if j==0:    
                    ax.set_xlabel('                     strength of A')
                    if len(vir_y)>0:    
                        ax.set_ylabel('strength of B')
                    else:
                        ax.set_yticks([])
                        # ax.set_ylabel('avirulent \n species B')
                fig.add_subplot(ax)
                
                plt.subplots_adjust(bottom=0.25, left=0.05, right=0.99, top=0.95)
                
                
    labels=['a)','b)','d)','e)']
    for i,l in zip([0,2,6,8],labels):
        ax=axes[i]
        axes[i].text(-0.3, 3.2, l, fontsize=8,transform=ax.transAxes, verticalalignment='top')
    labels=['c)','f)','g)','h)','i)'] 
    for i,l in zip([4,10,12,14],labels):
        ax=axes[i]
        axes[i].text(-0.3, 1.4, l, fontsize=8,transform=ax.transAxes, verticalalignment='top')
            
    # remove y ticks on species B            
    for k in range(1,len(axes),2):
        axes[k].set_yticks([])
    
    # replace ytick labels with L, M, H
    for k in range(len(axes)):
        if len(axes[k].get_ylabel())>0:
            if axes[k].get_ylabel()[0]=='s':
                axis_ticks_y = ['$H$', '$M$', '$L$']
                axes[k].set_yticks(np.arange(len(axis_ticks_y)), labels=axis_ticks_y)
            else:
                axes[k].set_yticklabels([])
    # x tick labels
    for k in range(len(axes)):
        axis_ticks_x = ['$L$', '$M$', '$H$']
        axes[k].set_xticks(np.arange(len(axis_ticks_x)), labels=axis_ticks_x)

    props = dict(boxstyle='round', facecolor='white', alpha=0.1,edgecolor='white')
    axes[3].text(-2.7,-1.2, 'Normalized AUC \n of each species', transform=ax.transAxes, verticalalignment='top', bbox=props)
    axes[3].text(3.3,-1.1, 'A\n \nB', transform=ax.transAxes, verticalalignment='top', bbox=props)

    saveloc = fig_path + "/AUC_v3s_"+str(v)
    fig.savefig(saveloc+".png", format="png",dpi=600)
     