#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
An agent based model of a polymicrobial community with virulence-mediated interspecies interactions
Author: Eske van Meegen and Cathi Herzberg, LACDR, Leiden University, The Netherlands

Parts of the code are based on Herzberg C. (2023) Vanhasseltlab/polymicro-abm: Polymicro-abm. v1.0.0. An agent
based modeling framework of an interactive, polymicrobial community during antimicrobial treatment. Zenodo. 
doi:10.5281/zenodo.10039816 which can be found on https://github.com/vanhasseltlab/PolyMicro-ABM 

File Description: This script creates Growthcurve plots of dual-species populations under treatment. Figures 8 and 9
(Filenames Fig14.png and Fig15.png)
"""
#%%
# %reset -f
# %matplotlib qt

import os
import pickle
import numpy as np
import pandas as pd    
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

#%% CALCULATE NEW MEAN - KEEPING ALL ITERATIONS UNTIL END OF RUNTIME - skip this and load database_filled.pkl if already created

#working directory: Main folder of repository "PolyMicro-ABM-Virulence"
path = os.getcwd() 
path_data = path + "/00_Complete_Database"
# load data
file = open(path_data +  "/database.pkl", "rb")
database = pickle.load(file)

#fill up short iterations
for key in database:
    for key2 in database[key]:
        for key3 in database[key][key2]:
            ls=[] # collect length of each iteration
            for key5 in database[key][key2][key3]['iterations']:
                l=len(database[key][key2][key3]['iterations'][key5])
                ls.append(l)
            maxl=max(ls) # maximal length of all iterations
            indexes=[i for i in range(len(ls)) if ls[i]!=maxl] #indexes of iterations that have not reached maximal lengths
            indexes_max=[i for i in range(len(ls)) if ls[i]==maxl]
            column_s=list(database[key][key2][key3]['iterations'][0].columns) # use column_names of first iteration for rows that are added
            for i in indexes: # go through dataframes that do not have maximal length
                for j in range(maxl-ls[i]): # for each missing timepoint
                    # add a row at the and of that dataframe
                    #determine timestep and Time from the longest iteration
                    timestep=database[key][key2][key3]['iterations'][indexes_max[0]].loc[ls[i]+j]['timestep']
                    Time=database[key][key2][key3]['iterations'][indexes_max[0]].loc[ls[i]+j]['Time']
                    #determine number of microbes etc. from the last recorded row of the dataframe that is being extended
                    Number_microbes=database[key][key2][key3]['iterations'][i].loc[ls[i]-1]['Number_microbes'] # last recorded timestep
                    Number_neutrophils=database[key][key2][key3]['iterations'][i].loc[ls[i]-1]['Number_neutrophils']
                    AB_fix_I=Number_neutrophils=database[key][key2][key3]['iterations'][i].loc[ls[i]-1]['AB_fix_I']
                    A=database[key][key2][key3]['iterations'][i].loc[ls[i]-1]['#A']
                    B=database[key][key2][key3]['iterations'][i].loc[ls[i]-1]['#B']
                    iteration=database[key][key2][key3]['iterations'][i].loc[ls[i]-1]['iteration']
                    #set to zero
                    Number_recruited=0
                    slope_A=0
                    slope_B=0
                    Bacterial_survival=0
                    #set the rest to nan 
                    n=len(database[key][key2][key3]['iterations'][i].columns)-12
                    database[key][key2][key3]['iterations'][i].loc[ls[i]+j] = [timestep,Time,Number_microbes,Number_neutrophils,Number_recruited,AB_fix_I,A,B,slope_A,
                                                                        slope_B,Bacterial_survival,iteration]+[np.nan]*n
                                                                                                                                
# calculate mean and add to mean df one level up as extra columns
for key in database:
    for key2 in database[key]:
        for key3 in database[key][key2]:
            #add iteration column to each dataframe
            for key5 in database[key][key2][key3]['iterations']:
                database[key][key2][key3]['iterations'][key5]['iteration']=key5
                l=len(database[key][key2][key3]['iterations'][key5])
                database[key][key2][key3]['iterations'][key5]['step']=list(range(l))
            #combine all iterations into one dataframe
            df=pd.concat(database[key][key2][key3]['iterations'],ignore_index=True,sort=False)
            #calculate mean across iterations
            col_names=['Number_microbes','#A','#B']
            df_m=df.groupby(['step'])[col_names].mean()
            df_m=df_m.rename({'Number_microbes':'M-Number_microbes','#A':'M-#A','#B':'M-#B'},axis=1)
            df_s=df.groupby(['step'])[col_names].std() 
            df_s=df_s.rename({'Number_microbes':'S-Number_microbes','#A':'S-#A','#B':'S-#B'},axis=1)
            #add to mean dataframe
            database[key][key2][key3]['mean']=pd.concat([database[key][key2][key3]['mean'],df_m,df_s],axis=1)
            
#save with newly calculated means
file = open(path_data + "/database_filled.pkl", "wb")
pickle.dump(database, file)
file.close()   
#%% load filled up database with newly calculated mean
file = open(path_data + "/database_filled.pkl", "rb")
database = pickle.load(file)

#%% plot

# define plot formatting settings
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['figure.figsize'] = [3.2,2]
mpl.rcParams['figure.dpi']= 600
mpl.rcParams['lines.linewidth']  = 1
mpl.rcParams['legend.fontsize'] = 7
mpl.rcParams['legend.title_fontsize'] = 8
mpl.rcParams['axes.labelsize']=8
mpl.rcParams['axes.titlesize']=8
mpl.rcParams['xtick.labelsize']=8 
mpl.rcParams['ytick.labelsize']=8
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['xtick.major.size'] = 1
mpl.rcParams['xtick.major.width'] = 0.5
mpl.rcParams['xtick.minor.size'] = 0.5
mpl.rcParams['xtick.minor.width'] = 0.3
mpl.rcParams['ytick.major.size'] = 1
mpl.rcParams['ytick.major.width'] = 0.5
mpl.rcParams['ytick.minor.size'] = 0.5
mpl.rcParams['ytick.minor.width'] = 0.3

# MAKE FIGURE FOLDER 
fig_folder = "/Figures/Growthcurves" 
fig_path = path + fig_folder
if not (os.path.exists(fig_path)): # if the folder doesn't exist yet
    os.mkdir(fig_path)
#%% FUNCTION TO SAVE LEGEND IN SEPERATE FIG
def save_legend(fig, l, saveloc, leg_title="Legend", leg_name='_Legend.png'):
    legend = fig.legend([x[0][0] for x in l],["{}".format(x[1]) for x in l], title=leg_title, bbox_to_anchor=(0.5, 0.5, 0, 0), ncol=1, loc='center', fontsize=8, title_fontsize=10)
    leg = legend.figure
    leg.canvas.draw()
    bbox = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(saveloc+leg_name, format="png", bbox_inches=bbox)    

#%% plot each scenario individually
# colors and linestyles
c = ['black', '#FFC107','#1E88E5','#D81B60','#009E73']
ls=['solid','dashed', 'dashdot']
runtime = 1440

for key in database:
    for key2 in database[key]:
        if len(database[key][key2])>1:
            fig, axes = plt.subplots(nrows=1,ncols=2,gridspec_kw={'wspace':0.03,'hspace':0.03}, sharex=True,sharey=True)
            c_i=0
            for key3 in database[key][key2]:
                x=database[key][key2][key3]['mean']['Time']
                for j in range(2):
                    ax=axes.flat[j]
                    species = key.split(' + ')[j][0]
                    y=database[key][key2][key3]['mean']['M-#'+ species]
                    y_s=database[key][key2][key3]['mean']['S-#'+ species]
                    ax.plot(x,y,color=c[c_i])
                    ax.fill_between(x,y-y_s,y+y_s,alpha=0.1,color=c[c_i])
                    vir_str = ' + '.join(['$v_{'+k+'}$' for k in key.split(' + ')[j][1::] if k!='_']) # creates v_ title 
                    if len(vir_str)>0:
                        vir_str += ' - '+ str(key2[j][0]) # add strength if species has virulence ability
                    ax.set_title('Species '+species+ '\n ' + vir_str) # generates title from key and key2
                    ax.set_xlabel('Time [min]')
                    ax.set_xlim([0,runtime])
                    ax.set_ylim([0,3500])
                    
                axes.flat[0].set_ylabel('Number of microbes')
                c_i+=1
            
            plt.subplots_adjust(bottom=0.2, left=0.2, right=0.99, top=0.85)
            saveloc = fig_path + "/Growthcurves_"+str(key)+"_"+str(key2)
            fig.savefig(saveloc+".png", format="png",dpi=600)


#%%
#legend, CORRECT THIS
fig, axes = plt.subplots(nrows=1,ncols=3,gridspec_kw={'wspace':0.05,'hspace':0.05}, sharex=True,sharey=True)
l2=[]
AB_concs = [0, 1, 1.15, 1.5, 4]
# make dummy lines for legend 
ax2 = ax.twinx() 
for i in range(len(AB_concs)):
    dummy = ax2.plot([], [], color=c[i], linestyle='solid',alpha=1) # for legend
    l2.append([dummy, str(AB_concs[i])]) # for legend
ax2.set_yticks([]) # removes ticks from second x-axis
save_legend(fig, l2, saveloc, leg_title="Concentration\n      (xMIC)", leg_name='_Legend.png')

#%% Plot growthcurves figures together for Figures Part 3 - Treatment
# define plot formatting settings
plt.style.use('default')
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['figure.dpi']= 600
mpl.rcParams['figure.figsize']=[3.2,3.2]
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

# colors and linestyles
c = ['black', '#FFC107','#1E88E5','#D81B60','#009E73']
ls=['solid','dashed', 'dashdot']
runtime = 1440


plot_keys=['Acs + B_', 'Ac + Bs', 'Acsr + B_', 'Acr + Bs', 'Acs + Br', 'Asr + Bc']
plot_keys2=[('High', 'None'),('High', 'High'),('High', 'None'),('High', 'High'),('High', 'High'),('High', 'High')]

fig = plt.figure()
axes=[]
outer = gridspec.GridSpec(3, 2, wspace=0.3, hspace=0.8)
i=0
for key,key2,i in zip(plot_keys,plot_keys2,range(6)):
    inner = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[i], wspace=0.1, hspace=0)
    
    for j in range(2):
        x_allconcs=[database[key][key2][key3]['mean']['Time'] for key3 in database[key][key2]]
        y_allconcs=[database[key][key2][key3]['mean']['M-#'+ key.split(' + ')[j][0]] for key3 in database[key][key2]]
        ys_allconcs=[database[key][key2][key3]['mean']['S-#'+ key.split(' + ')[j][0]] for key3 in database[key][key2]]
        concs=[key3 for key3 in database[key][key2]]
        ax = plt.Subplot(fig, inner[j])
        axes=axes+[ax]
        c_i=0
        for x,y,y_s,conc in zip(x_allconcs,y_allconcs,ys_allconcs,concs):
            ax.plot(x,y,color=c[c_i])
            ax.fill_between(x,y-y_s,y+y_s,alpha=0.1,color=c[c_i])
            c_i+=1
        vir_str = ' + '.join(['$v_{'+k+'}$' for k in key.split(' + ')[j][1::] if k!='_']) # creates v_ title 
        species = key.split(' + ')[j][0]
        ax.set_title(species+ '\n ' + vir_str) # generates title from key and key2
        
        ax.set_xlim([0,runtime])
        ax.set_ylim([0,3500])
        ax.grid(which='major')
        fig.add_subplot(ax)

axes[8].set_xlabel('               Time [min]')
axes[10].set_xlabel('               Time [min]')

axes[0].set_ylabel('Number of \n microbes')
axes[4].set_ylabel('Number of \n microbes')
axes[8].set_ylabel('Number of \n microbes')
     
axes[1].set_yticklabels([])
axes[3].set_yticklabels([])
axes[5].set_yticklabels([])
axes[7].set_yticklabels([])
axes[9].set_yticklabels([])
axes[11].set_yticklabels([])      
axes[2].set_yticklabels([])  
axes[6].set_yticklabels([])  
axes[10].set_yticklabels([]) 

for i in range(8):
    axes[i].set_xticklabels([]) 

y=1.5
labels=['a)','c)','e)']
for i,l in zip([0,4,8],labels):
    ax=axes[i]
    axes[i].text(-0.35, y, l, fontsize=8,transform=ax.transAxes, verticalalignment='top')

labels=['b)','d)','f)']
for i,l in zip([2,6,10],labels):
    ax=axes[i]
    axes[i].text(-0.35, y, l, fontsize=8,transform=ax.transAxes, verticalalignment='top')

#add legend
legend_elements = [
    Line2D([0], [0], color=c[0], label=concs[0]),
    Line2D([0], [0], color=c[1], label=concs[1]), 
    Line2D([0], [0], color=c[2], label=concs[2]), 
    Line2D([0], [0], color=c[3], label=concs[3]), 
    Line2D([0], [0], color=c[4], label=concs[4])]
ax=axes[10]
ax.legend(handles=legend_elements,title='Drug concentration',ncol=5,loc='lower right', bbox_to_anchor=(1.5, -1.25))
plt.subplots_adjust(bottom=0.19, left=0.2, right=0.98, top=0.9)

saveloc = fig_path + "/Fig14"
fig.savefig(saveloc+".png", format="png",dpi=600)

#%%
mpl.rcParams['figure.figsize']=[3.2,2.4]

plot_keys=['As + B_', 'Asr + B_', 'As + Br', 'Asr + Br']
plot_keys2=[('High', 'None'),('High', 'None'),('High', 'High'),('High', 'High')]

fig = plt.figure()
axes=[]
outer = gridspec.GridSpec(2, 2, wspace=0.3, hspace=0.8)
i=0
for key,key2,i in zip(plot_keys,plot_keys2,range(4)):
    inner = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[i], wspace=0.1, hspace=0)
    
    for j in range(2):
        x_allconcs=[database[key][key2][key3]['mean']['Time'] for key3 in database[key][key2]]
        y_allconcs=[database[key][key2][key3]['mean']['M-#'+ key.split(' + ')[j][0]] for key3 in database[key][key2]]
        ys_allconcs=[database[key][key2][key3]['mean']['S-#'+ key.split(' + ')[j][0]] for key3 in database[key][key2]]
        concs=[key3 for key3 in database[key][key2]]
        ax = plt.Subplot(fig, inner[j])
        axes=axes+[ax]
        c_i=0
        for x,y,y_s,conc in zip(x_allconcs,y_allconcs,ys_allconcs,concs):
            ax.plot(x,y,color=c[c_i])
            ax.fill_between(x,y-y_s,y+y_s,alpha=0.1,color=c[c_i])
            c_i+=1
        vir_str = ' + '.join(['$v_{'+k+'}$' for k in key.split(' + ')[j][1::] if k!='_']) # creates v_ title 
        species = key.split(' + ')[j][0]
        ax.set_title(species+ '\n ' + vir_str) # generates title from key and key2
        
        ax.set_xlim([0,runtime])
        ax.set_ylim([0,3500])
        ax.grid(which='major')
        fig.add_subplot(ax)

axes[4].set_xlabel('               Time [min]')
axes[6].set_xlabel('               Time [min]')

axes[0].set_ylabel('Number of \n microbes')
axes[4].set_ylabel('Number of \n microbes')
     
axes[1].set_yticklabels([])
axes[3].set_yticklabels([])
axes[5].set_yticklabels([])
axes[7].set_yticklabels([])
axes[2].set_yticklabels([])  
axes[6].set_yticklabels([])  

for i in range(4):
    axes[i].set_xticklabels([]) 

y=1.5
labels=['a)','b)','c)','d)']
for i,l in zip([0,2,4,6],labels):
    ax=axes[i]
    axes[i].text(-0.35, y, l, fontsize=8,transform=ax.transAxes, verticalalignment='top')


#add legend
legend_elements = [
    Line2D([0], [0], color=c[0], label=concs[0]),
    Line2D([0], [0], color=c[1], label=concs[1]), 
    Line2D([0], [0], color=c[2], label=concs[2]), 
    Line2D([0], [0], color=c[3], label=concs[3]), 
    Line2D([0], [0], color=c[4], label=concs[4])]
ax=axes[6]
ax.legend(handles=legend_elements,title='Drug concentration',ncol=5,loc='lower right', bbox_to_anchor=(1.5, -1.25))
plt.subplots_adjust(bottom=0.27, left=0.2, right=0.98, top=0.85)

saveloc = fig_path + "/Fig15"
fig.savefig(saveloc+".png", format="png",dpi=600)
