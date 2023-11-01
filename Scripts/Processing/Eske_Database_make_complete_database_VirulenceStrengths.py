# -*- coding: utf-8 -*-
"""
An agent based model of a polymicrobial community with virulence-mediated interspecies interactions
Author: Eske van Meegen and Cathi Herzberg, LACDR, Leiden University, The Netherlands

Parts of the code are based on Herzberg C. (2023) Vanhasseltlab/polymicro-abm: Polymicro-abm. v1.0.0. An agent
based modeling framework of an interactive, polymicrobial community during antimicrobial treatment. Zenodo. 
doi:10.5281/zenodo.10039816 which can be found on https://github.com/vanhasseltlab/PolyMicro-ABM 

File Description: This script creates a database which includes simulation data from all simulations with a range of virulence strengths
"""
#%%
import os
import pickle
import numpy as np
import itertools
import time
import pandas as pd
import math

#%% OPEN 03_VIRULENCE_STRENGTHS_A1+B0
path_core = os.getcwd() + "/ownCloud2/cathi/Polymicro-not-shared/Virulence paper with Eske/Data/" 

path = path_core + "03_Virulence_Strengths_A1+B0/220612_vf_cellbound"

file = open(path + "/db_m_part.pkl", "rb")
db_m_part = pickle.load(file)

file = open(path + "/dict_model.pkl", "rb")
dict_model = pickle.load(file)

p=1 #position of virulence strengths in db_m_part key

#%% MAKE DATABASE VC DICT for dual species 

# insert virtype A + virtype B
database2c = {}
virtype='Ac + B_'
conc=0
database2c[virtype] = {}

# insert strength A + strength B per virtypes
virstrengths=[(str(y),'None') for y in list(dict.fromkeys([x[p] for x in db_m_part.keys()]))]
for k in virstrengths:
    database2c[virtype][k] = {}
    database2c[virtype][k][conc]={} # insert AB concentrations per strengths
    dicts = ['settings', 'mean', 'iterations']
    for i in dicts:
        database2c[virtype][k][conc][i] = {} # insert settings, mean and iterations per ab_conc
            
#%% ADD MEANS TO DATABASE
# from dict_model_pro into database 'mean'

for key in dict_model:
    virstrengths=(str(key[p]),'None')
    database2c[virtype][virstrengths][conc]['mean'] = dict_model[key]

#%% ADD ITERATIONS TO DATABASE
# from dict_model_pro into database 'iterations'
for key in db_m_part:
    virstrengths=(str(key[p]),'None')
    it=key[-1]
    database2c[virtype][virstrengths][conc]['iterations'][it]=db_m_part[key]
#%%













#%% OPEN 03_VIRULENCE_STRENGTHS_A1+B0
path_core = os.getcwd() + "/ownCloud2/cathi/Polymicro-not-shared/Virulence paper with Eske/Data/" 

path = path_core + "03_Virulence_Strengths_A1+B0/220612_vf_secreted"

file = open(path + "/db_m_part.pkl", "rb")
db_m_part = pickle.load(file)

file = open(path + "/dict_model.pkl", "rb")
dict_model = pickle.load(file)

p=2 #position of virulence strengths in db_m_part key
#%% MAKE DATABASE VS DICT for dual species 

# insert virtype A + virtype B
database2s = {}
virtype='As + B_'
conc=0
database2s[virtype] = {}

# insert strength A + strength B per virtypes
virstrengths=[(str(y),'None') for y in list(dict.fromkeys([x[p] for x in db_m_part.keys()]))]
for k in virstrengths:
    database2s[virtype][k] = {}
    database2s[virtype][k][conc]={} # insert AB concentrations per strengths
    dicts = ['settings', 'mean', 'iterations']
    for i in dicts:
        database2s[virtype][k][conc][i] = {} # insert settings, mean and iterations per ab_conc
            
#%% ADD MEANS TO DATABASE
# from dict_model_pro into database 'mean'

for key in dict_model:
    virstrengths=(str(key[p]),'None')
    database2s[virtype][virstrengths][conc]['mean'] = dict_model[key]

#%% ADD ITERATIONS TO DATABASE
# from dict_model_pro into database 'iterations'
for key in db_m_part:
    virstrengths=(str(key[p]),'None')
    it=key[-1]
    database2s[virtype][virstrengths][conc]['iterations'][it]=db_m_part[key]    
 
#%%













#%% OPEN 03_VIRULENCE_STRENGTHS_A1+B0
path_core = os.getcwd() + "/ownCloud2/cathi/Polymicro-not-shared/Virulence paper with Eske/Data/" 

path = path_core + "03_Virulence_Strengths_A1+B0/220612_vf_recruit"

file = open(path + "/db_m_part.pkl", "rb")
db_m_part = pickle.load(file)

file = open(path + "/dict_model.pkl", "rb")
dict_model = pickle.load(file)

p=3 #position of virulence strengths in db_m_part key

#%% MAKE DATABASE VS DICT for dual species 

# insert virtype A + virtype B
database2r = {}
virtype='Ar + B_'
conc=0
database2r[virtype] = {}

# insert strength A + strength B per virtypes
virstrengths=[(str(y),'None') for y in list(dict.fromkeys([x[p] for x in db_m_part.keys()]))]
for k in virstrengths:
    database2r[virtype][k] = {}
    database2r[virtype][k][conc]={} # insert AB concentrations per strengths
    dicts = ['settings', 'mean', 'iterations']
    for i in dicts:
        database2r[virtype][k][conc][i] = {} # insert settings, mean and iterations per ab_conc
            
#%% ADD MEANS TO DATABASE
# from dict_model_pro into database 'mean'

for key in dict_model:
    virstrengths=(str(key[p]),'None')
    database2r[virtype][virstrengths][conc]['mean'] = dict_model[key]

#%% ADD ITERATIONS TO DATABASE
# from dict_model_pro into database 'iterations'
for key in db_m_part:
    virstrengths=(str(key[p]),'None')
    it=key[-1]
    database2r[virtype][virstrengths][conc]['iterations'][it]=db_m_part[key]  
#%%
#
#
#
#
#
#
#
#
#
#
#

#%%













#%% OPEN 02_VIRULENCE_STRENGTHS_A1
path_core = os.getcwd() + "/ownCloud2/cathi/Polymicro-not-shared/Virulence paper with Eske/Data/" 

path = path_core + "02_Virulence_Strengths_A1/220527_vf_cellbound"

file = open(path + "/db_m_part.pkl", "rb")
db_m_part = pickle.load(file)

file = open(path + "/dict_model.pkl", "rb")
dict_model = pickle.load(file)

p=1 #position of virulence strengths in db_m_part key

#%% MAKE DATABASE VS DICT for dual species 

# insert virtype A + virtype B
database1c = {}
virtype='Ac'
conc=0
database1c[virtype] = {}

# insert strength A + strength B per virtypes
virstrengths=[(str(y),'None') for y in list(dict.fromkeys([x[p] for x in db_m_part.keys()]))]
for k in virstrengths:
    database1c[virtype][k] = {}
    database1c[virtype][k][conc]={} # insert AB concentrations per strengths
    dicts = ['settings', 'mean', 'iterations']
    for i in dicts:
        database1c[virtype][k][conc][i] = {} # insert settings, mean and iterations per ab_conc
            
#%% ADD MEANS TO DATABASE
# from dict_model_pro into database 'mean'

for key in dict_model:
    virstrengths=(str(key[p]),'None')
    database1c[virtype][virstrengths][conc]['mean'] = dict_model[key]

#%% ADD ITERATIONS TO DATABASE
# from dict_model_pro into database 'iterations'
for key in db_m_part:
    virstrengths=(str(key[p]),'None')
    it=key[-1]
    database1c[virtype][virstrengths][conc]['iterations'][it]=db_m_part[key]    
#%% add second set of iterations in 02_VIRULENCE_STRENGTHS_A1
path_core = os.getcwd() + "/ownCloud2/cathi/Polymicro-not-shared/Virulence paper with Eske/Data/" 

path = path_core + "02_Virulence_Strengths_A1/220622_vf_cellbound"

file = open(path + "/db_m_part.pkl", "rb")
db_m_part = pickle.load(file)

file = open(path + "/dict_model.pkl", "rb")
dict_model = pickle.load(file)

p=1 #position of virulence strengths in db_m_part key

#%% ADD ITERATIONS TO DATABASE
# from dict_model_pro into database 'iterations'
any_key=list(database1c[virtype].keys())[0]
n_its=len(database1c[virtype][any_key][conc]['iterations']) # number of its already in database
for key in db_m_part:
    virstrengths=(str(key[p]),'None')
    it=key[-1]+n_its
    database1c[virtype][virstrengths][conc]['iterations'][it]=db_m_part[key] 
            
#%% ADD MEANS TO DATABASE
# from dict_model_pro into database 'mean'

# for key in dict_model:
#     virstrengths=(str(key[p]),'None')
#     database1r[virtype][virstrengths][conc]['mean'] = dict_model[key]
#%%













#%% OPEN 02_VIRULENCE_STRENGTHS_A1
path_core = os.getcwd() + "/ownCloud2/cathi/Polymicro-not-shared/Virulence paper with Eske/Data/" 

path = path_core + "02_Virulence_Strengths_A1/220527_vf_secreted"

file = open(path + "/db_m_part.pkl", "rb")
db_m_part = pickle.load(file)

file = open(path + "/dict_model.pkl", "rb")
dict_model = pickle.load(file)

p=2 #position of virulence strengths in db_m_part key

#%% MAKE DATABASE VS DICT for dual species 

# insert virtype A + virtype B
database1s = {}
virtype='As'
conc=0
database1s[virtype] = {}

# insert strength A + strength B per virtypes
virstrengths=[(str(y),'None') for y in list(dict.fromkeys([x[p] for x in db_m_part.keys()]))]
for k in virstrengths:
    database1s[virtype][k] = {}
    database1s[virtype][k][conc]={} # insert AB concentrations per strengths
    dicts = ['settings', 'mean', 'iterations']
    for i in dicts:
        database1s[virtype][k][conc][i] = {} # insert settings, mean and iterations per ab_conc
            
#%% ADD MEANS TO DATABASE
# from dict_model_pro into database 'mean'

# for key in dict_model:
#     virstrengths=(str(key[p]),'None')
#     database1s[virtype][virstrengths][conc]['mean'] = dict_model[key]

#%% ADD ITERATIONS TO DATABASE
# from dict_model_pro into database 'iterations'
for key in db_m_part:
    virstrengths=(str(key[p]),'None')
    it=key[-1]
    database1s[virtype][virstrengths][conc]['iterations'][it]=db_m_part[key]    
#%% add second set of iterations in 02_VIRULENCE_STRENGTHS_A1
path_core = os.getcwd() + "/ownCloud2/cathi/Polymicro-not-shared/Virulence paper with Eske/Data/" 

path = path_core + "02_Virulence_Strengths_A1/220622_vf_secreted"

file = open(path + "/db_m_part.pkl", "rb")
db_m_part = pickle.load(file)

file = open(path + "/dict_model.pkl", "rb")
dict_model = pickle.load(file)

p=2 #position of virulence strengths in db_m_part key

#%% ADD ITERATIONS TO DATABASE
# from dict_model_pro into database 'iterations'
any_key=list(database1s[virtype].keys())[0]
n_its=len(database1s[virtype][any_key][conc]['iterations']) # number of its already in database
for key in db_m_part:
    virstrengths=(str(key[p]),'None')
    if virstrengths in database1s[virtype].keys():
        it=key[-1]+n_its
        database1s[virtype][virstrengths][conc]['iterations'][it]=db_m_part[key]
    else:
        database1s[virtype][virstrengths] = {}
        database1s[virtype][virstrengths][conc]={} # insert AB concentrations per strengths
        dicts = ['settings', 'mean', 'iterations']
        for i in dicts:
            database1s[virtype][virstrengths][conc][i] = {} # insert settings, mean and iterations per ab_conc
        it=key[-1]
        database1s[virtype][virstrengths][conc]['iterations'][it]=db_m_part[key]   
        
            
#%% ADD MEANS TO DATABASE
# from dict_model_pro into database 'mean'

for key in dict_model:
    virstrengths=(str(key[p]),'None')
    database1s[virtype][virstrengths][conc]['mean'] = dict_model[key]

#%%













#%% OPEN 02_VIRULENCE_STRENGTHS_A1
path_core = os.getcwd() + "/ownCloud2/cathi/Polymicro-not-shared/Virulence paper with Eske/Data/" 

path = path_core + "02_Virulence_Strengths_A1/220607_vf_recruit"

file = open(path + "/db_m_part.pkl", "rb")
db_m_part = pickle.load(file)

file = open(path + "/dict_model.pkl", "rb")
dict_model = pickle.load(file)

p=3 #position of virulence strengths in db_m_part key

#%% MAKE DATABASE VS DICT for dual species 

# insert virtype A + virtype B
database1r = {}
virtype='Ar'
conc=0
database1r[virtype] = {}

# insert strength A + strength B per virtypes
virstrengths=[(str(y),'None') for y in list(dict.fromkeys([x[p] for x in db_m_part.keys()]))]
for k in virstrengths:
    database1r[virtype][k] = {}
    database1r[virtype][k][conc]={} # insert AB concentrations per strengths
    dicts = ['settings', 'mean', 'iterations']
    for i in dicts:
        database1r[virtype][k][conc][i] = {} # insert settings, mean and iterations per ab_conc
            
#%% ADD MEANS TO DATABASE
# from dict_model_pro into database 'mean'

# for key in dict_model:
#     virstrengths=(str(key[p]),'None')
#     database1r[virtype][virstrengths][conc]['mean'] = dict_model[key]

#%% ADD ITERATIONS TO DATABASE
# from dict_model_pro into database 'iterations'
for key in db_m_part:
    virstrengths=(str(key[p]),'None')
    it=key[-1]
    database1r[virtype][virstrengths][conc]['iterations'][it]=db_m_part[key]    
#%% add second set of iterations in 02_VIRULENCE_STRENGTHS_A1
path_core = os.getcwd() + "/ownCloud2/cathi/Polymicro-not-shared/Virulence paper with Eske/Data/" 

path = path_core + "02_Virulence_Strengths_A1/220622_vf_recruit"

file = open(path + "/db_m_part.pkl", "rb")
db_m_part = pickle.load(file)

file = open(path + "/dict_model.pkl", "rb")
dict_model = pickle.load(file)

p=3 #position of virulence strengths in db_m_part key

#%% ADD ITERATIONS TO DATABASE
# from dict_model_pro into database 'iterations'
any_key=list(database1r[virtype].keys())[0]
n_its=len(database1r[virtype][any_key][conc]['iterations']) # number of its already in database
for key in db_m_part:
    if key[p]==0 or key[p]==1:
        virstrengths=(str(int(key[p])),'None')
    else:
        virstrengths=(str(key[p]),'None')
    it=key[-1]+n_its
    database1r[virtype][virstrengths][conc]['iterations'][it]=db_m_part[key] 
            
#%% ADD MEANS TO DATABASE
# from dict_model_pro into database 'mean'

for key in dict_model:
    if key[p]==0 or key[p]==1:
        virstrengths=(str(int(key[p])),'None')
    else:
        virstrengths=(str(key[p]),'None')
    
    database1r[virtype][virstrengths][conc]['mean'] = dict_model[key]
#%%
#
#
#
#
#
#
#
#
#
#
#%% MERDE DATABASE93, DATABASE94 AND DATABASE95 INTO COMPLETE DATABASE

database = database2c | database2s | database2r | database1c | database1s | database1r

#%% SAVE COMPLETE DATABASE

# os.chdir(r'D:\Research Project 1')
path = path_core + "001_Complete_Database_Virulence_Strengths"
file = open(path + "/database.pkl", "wb")
pickle.dump(database, file)
file.close() 
#%%
import os
import pickle
import numpy as np
import itertools
import time
import pandas as pd
import math

#%% OPEN COMPLETE DATABASE
path_core = os.getcwd() + "/ownCloud2/cathi/Polymicro-not-shared/Virulence paper with Eske/Data/" 
path = path_core + "001_Complete_Database_Virulence_Strengths"
file = open(path + "/database.pkl", "rb")
database = pickle.load(file)
#%%

"""
Add settings 
"""

#%% DEF FUNCTIONS
        
def open_all_settings(folder_path, folder_name):
    
    # change to correct folder
    path = path_core + folder_path +'/' + folder_name
    
    # open files
    file = open(path + "/species.pkl", "rb")
    all_species,all_species_all_virulence,ab_conc_range = pickle.load(file)
    file.close()

    file = open(path + "/neutrophils.pkl", "rb")
    neutrophils = pickle.load(file)
    file.close()

    file = open(path + "/settings.pkl", "rb")
    microbe_capacity,n_init,neutrophil_capacity,neutrophil_n_init,its,runtime,gridrange,stepsize,data_collection,variable_params,fixed_params = pickle.load(file)
    file.close()
    
    # collect all settings
    max_recruit_rate = neutrophils[0][0]
    halflife = neutrophils[0][1]
    phag_rate = neutrophils[0][2]

    drug_Gmin = fixed_params['drugs'][0]
    drug_k = fixed_params['drugs'][1]
    drug_type = fixed_params['drugs'][2]

    all_settings = {
                    # general settings
                    'its': its,
                    'runtime': runtime,
                    'gridrange': gridrange,
                    'stepsize': stepsize,
                    'data_collection': data_collection,
                    
                    # microbes
                    'microbe_capacity': microbe_capacity,
                    'n_init': n_init,
                    'all_species': all_species,
                    #'all_species_all_virulence': all_species_all_virulence,
                    
                    # drug
                    'drug_Gmin': drug_Gmin,
                    'drug_k': drug_k,
                    'drug_type': drug_type,
                    #'ab_conc_range': ab_conc_range,
                    
                    # neutrophils
                    'neutrophil_capacity': neutrophil_capacity,
                    'neutrophil_n_init': neutrophil_n_init,                
                    'max_recruit_rate': max_recruit_rate,
                    'halflife': halflife,
                    'phag_rate': phag_rate
                    }
    
    return all_settings

def add_settings(database, all_settings, virtype, abconc, virstrength=None):
    if abconc == 0:
        if virstrength == None:
            for key in database[virtype]:
                database[virtype][key][abconc]['settings'] = all_settings
        else:
            database[virtype][virstrength][abconc]['settings'] = all_settings  
            
    if abconc == [5, 5.75, 7, 20]:
        for conc in abconc:          
            if virstrength == None:
                for key in database[virtype]:
                    database[virtype][key][conc]['settings'] = all_settings
            else:
                database[virtype][virstrength][conc]['settings'] = all_settings 
            
        
#%%

abconcs = [0]

# [[folder_path, folder_name, all_settings, virtype, abconc, virstrength]]

datas = [["03_Virulence_Strengths_A1+B0", "220612_vf_cellbound", 'Ac + B_', 0, None],
         ["03_Virulence_Strengths_A1+B0", "220612_vf_secreted", 'As + B_', 0, None],
         ["03_Virulence_Strengths_A1+B0", "220612_vf_recruit", 'Ar + B_', 0, None],
         
         ["02_Virulence_Strengths_A1", "220622_vf_cellbound", 'Ac', 0, None],
         ["02_Virulence_Strengths_A1", "220622_vf_secreted", 'As', 0, None],
         ["02_Virulence_Strengths_A1", "220622_vf_recruit", 'Ar', 0, None],
        ]

for data in datas:
    all_settings = open_all_settings(data[0], data[1])
    add_settings(database, all_settings, data[2], data[3], data[4])
    
#%% check if all simulations have their settings added

for virtype in database:
    for virstrength in database[virtype]:
        for abconc in database[virtype][virstrength]:
            booleans = []
            #print(len(database[virtype][virstrength][abconc]['mean']))
            if len(database[virtype][virstrength][abconc]['mean']) == 0:
                booleans.append(False)
            if len(database[virtype][virstrength][abconc]['iterations']) == 0:
                booleans.append(False)
            if len(database[virtype][virstrength][abconc]['settings']) == 0:
                booleans.append(False)
            
            if len(booleans) == 0: # if dict is not empty
                exit
            elif len(booleans) == 3:
                print()
            else:
                print('ERROR')
                print(len(booleans))
                print(virtype, virstrength, abconc)
                

#%% SAVE COMPLETE DATABASE

# os.chdir(r'D:\Research Project 1')
path = path_core + "001_Complete_Database_Virulence_Strengths"
file = open(path + "/database.pkl", "wb")
pickle.dump(database, file)
file.close() 
#%%
import os
import pickle
import numpy as np
import itertools
import time
import pandas as pd
import math

#%% OPEN COMPLETE DATABASE
path_core = os.getcwd() + "/ownCloud2/cathi/Polymicro-not-shared/Virulence paper with Eske/Data/" 
path = path_core + "001_Complete_Database_Virulence_Strengths"
file = open(path + "/database.pkl", "rb")
database = pickle.load(file)

#%% CALCULATE NEW MEAN - KEEPING ALL ITERATIONS UNTIL END OF RUNTIME - skip this and load database_filled.pkl if already created


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
                    if 'B' in key:
                        B=database[key][key2][key3]['iterations'][i].loc[ls[i]-1]['#B']
                    else:
                        B=None
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
            if 'B' in key:
                col_names=['Number_microbes','#A','#B']
            else:
                col_names=['Number_microbes','#A']
            df_m=df.groupby(['step'])[col_names].mean()
            if 'B' in key:
               df_m=df_m.rename({'Number_microbes':'M-Number_microbes','#A':'M-#A','#B':'M-#B'},axis=1)
            else:
                df_m=df_m.rename({'Number_microbes':'M-Number_microbes','#A':'M-#A'},axis=1)
            df_s=df.groupby(['step'])[col_names].std() 
            if 'B' in key:
                df_s=df_s.rename({'Number_microbes':'S-Number_microbes','#A':'S-#A','#B':'S-#B'},axis=1)
            else:
                df_s=df_s.rename({'Number_microbes':'S-Number_microbes','#A':'S-#A'},axis=1)
            #add to mean dataframe
            database[key][key2][key3]['mean']=pd.concat([database[key][key2][key3]['mean'],df_m,df_s],axis=1)
            
#save with newly calculated means
file = open(path + "/database_filled.pkl", "wb")
pickle.dump(database, file)
file.close()   
#%% load filled up database with newly calculated mean
file = open(path + "/database_filled.pkl", "rb")
database = pickle.load(file)
      
