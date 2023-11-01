# -*- coding: utf-8 -*-
"""
An agent based model of a polymicrobial community with virulence-mediated interspecies interactions
Author: Eske van Meegen and Cathi Herzberg, LACDR, Leiden University, The Netherlands

Parts of the code are based on Herzberg C. (2023) Vanhasseltlab/polymicro-abm: Polymicro-abm. v1.0.0. An agent
based modeling framework of an interactive, polymicrobial community during antimicrobial treatment. Zenodo. 
doi:10.5281/zenodo.10039816 which can be found on https://github.com/vanhasseltlab/PolyMicro-ABM 

File Description: This script creates a database which includes simulation data from all dualspecies simulations with numeric virulence input (L,M,H)
"""
#%%

import os
import pickle
import numpy as np
import itertools
import time
import pandas as pd
import math

#%% OPEN 93_VIRULENCE

# A1/2/3 + B0 
os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\93_Virulence_A123+B0\\"

file = open(path + "/db_m_part_pro.pkl", "rb")
db_m_part_pro = pickle.load(file)

file = open(path + "/dict_model_pro.pkl", "rb")
dict_model_pro = pickle.load(file)

#%% MAKE DATABASE93 DICT

# insert virtype A + virtype B
database93 = {}
for virtype in db_m_part_pro:
    database93[virtype] = {}

# insert strength A + strength B per virtypes
for virtype in db_m_part_pro:
    for key in db_m_part_pro[virtype]:
        virstrengths = []
        if key[:2] not in virstrengths:
            virstrengths.append(key[:2])
        for k in virstrengths:
            database93[virtype][k] = {}
 
ab_concs = [0, 5, 5.75, 7.0, 20]
dicts = ['settings', 'mean', 'iterations']

for vir_type in db_m_part_pro:
    for key in db_m_part_pro[vir_type]:
        strength = key[:2]
        
        # insert AB concentrations per strengths 
        for conc in ab_concs:
            database93[vir_type][strength][conc] = {}
            
            # insert settings, mean and iterations per ab_conc
            for i in dicts:
                database93[vir_type][strength][conc][i] = {}
            
#%% ADD MEANS TO DATABASE
# from dict_model_pro into database 'mean'

for virtype in dict_model_pro:
    for virstrength_abconc in dict_model_pro[virtype]:
        virstrength = virstrength_abconc[:-1]
        abconc = virstrength_abconc[-1]
        
        database93[virtype][virstrength][abconc]['mean'] = dict_model_pro[virtype][virstrength_abconc]
        
#%% ADD ITERATIONS TO DATABASE
# from dict_model_pro into database 'iterations'

for virtype in db_m_part_pro:
    for virstrength_abconc_it in db_m_part_pro[virtype]:
        virstrength = virstrength_abconc_it[:2]
        abconc = virstrength_abconc_it[2]
        it = virstrength_abconc_it[-1]
        
        database93[virtype][virstrength][abconc]['iterations'][it] = db_m_part_pro[virtype][virstrength_abconc_it]

#%% SAFE FILE - DATABASE 93

os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\00_Complete_Database\\"

file = open(path + "/database93.pkl", "wb")
pickle.dump(database93, file)
file.close()  

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
%reset -f
%matplotlib qt

import os
import pickle
import numpy as np
import itertools
import time
import pandas as pd
import math

#%% OPEN 94_VIRULENCE

# A1/2/3 + B0 
os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\94_Virulence_A1+B1\\"

file = open(path + "/db_m_part_pro.pkl", "rb")
db_m_part_pro = pickle.load(file)

file = open(path + "/dict_model_pro.pkl", "rb")
dict_model_pro = pickle.load(file)

#%% MAKE DATABASE94 DICT

# insert virtype A + virtype B
database94 = {}
for virtype in db_m_part_pro:
    database94[virtype] = {}

# insert strength A + strength B per virtypes
for virtype in db_m_part_pro:
    for key in db_m_part_pro[virtype]:
        virstrengths = []
        if key[:2] not in virstrengths:
            virstrengths.append(key[:2])
        for k in virstrengths:
            database94[virtype][k] = {}
 
ab_concs = [0, 5, 5.75, 7.0, 20]
dicts = ['settings', 'mean', 'iterations']

for vir_type in db_m_part_pro:
    for key in db_m_part_pro[vir_type]:
        strength = key[:2]
        
        # insert AB concentrations per strengths 
        for conc in ab_concs:
            database94[vir_type][strength][conc] = {}
            
            # insert settings, mean and iterations per ab_conc
            for i in dicts:
                database94[vir_type][strength][conc][i] = {}
            
#%% ADD MEANS TO DATABASE
# from dict_model_pro into database 'mean'

for virtype in dict_model_pro:
    for virstrength_abconc in dict_model_pro[virtype]:
        virstrength = virstrength_abconc[:-1]
        abconc = virstrength_abconc[-1]
        
        database94[virtype][virstrength][abconc]['mean'] = dict_model_pro[virtype][virstrength_abconc]
        
#%% ADD ITERATIONS TO DATABASE
# from dict_model_pro into database 'iterations'

for virtype in db_m_part_pro:
    for virstrength_abconc_it in db_m_part_pro[virtype]:
        virstrength = virstrength_abconc_it[:2]
        abconc = virstrength_abconc_it[2]
        it = virstrength_abconc_it[-1]
        
        database94[virtype][virstrength][abconc]['iterations'][it] = db_m_part_pro[virtype][virstrength_abconc_it]

#%% SAFE FILE - DATABASE 94

os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\00_Complete_Database\\"

file = open(path + "/database94.pkl", "wb")
pickle.dump(database94, file)
file.close()  

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
%reset -f
%matplotlib qt

import os
import pickle
import numpy as np
import itertools
import time
import pandas as pd
import math

#%% OPEN 95_VIRULENCE

# A1/2/3 + B0 
os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\95_Virulence_A2+B1\\"

file = open(path + "/db_m_part_pro.pkl", "rb")
db_m_part_pro = pickle.load(file)

file = open(path + "/dict_model_pro.pkl", "rb")
dict_model_pro = pickle.load(file)

#%% MAKE DATABASE95 DICT

# insert virtype A + virtype B
database95 = {}
for virtype in db_m_part_pro:
    database95[virtype] = {}

# insert strength A + strength B per virtypes
for virtype in db_m_part_pro:
    for key in db_m_part_pro[virtype]:
        virstrengths = []
        if key[:2] not in virstrengths:
            virstrengths.append(key[:2])
        for k in virstrengths:
            database95[virtype][k] = {}
 
ab_concs = [0, 5, 5.75, 7.0, 20]
dicts = ['settings', 'mean', 'iterations']

for vir_type in db_m_part_pro:
    for key in db_m_part_pro[vir_type]:
        strength = key[:2]
        
        # insert AB concentrations per strengths 
        for conc in ab_concs:
            database95[vir_type][strength][conc] = {}
            
            # insert settings, mean and iterations per ab_conc
            for i in dicts:
                database95[vir_type][strength][conc][i] = {}
            
#%% ADD MEANS TO DATABASE
# from dict_model_pro into database 'mean'

for virtype in dict_model_pro:
    for virstrength_abconc in dict_model_pro[virtype]:
        virstrength = virstrength_abconc[:-1]
        abconc = virstrength_abconc[-1]
        
        database95[virtype][virstrength][abconc]['mean'] = dict_model_pro[virtype][virstrength_abconc]
        
#%% ADD ITERATIONS TO DATABASE
# from dict_model_pro into database 'iterations'

for virtype in db_m_part_pro:
    for virstrength_abconc_it in db_m_part_pro[virtype]:
        virstrength = virstrength_abconc_it[:2]
        abconc = virstrength_abconc_it[2]
        it = virstrength_abconc_it[-1]
        
        database95[virtype][virstrength][abconc]['iterations'][it] = db_m_part_pro[virtype][virstrength_abconc_it]

#%% SAFE FILE - DATABASE 95

os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\00_Complete_Database\\"

file = open(path + "/database95.pkl", "wb")
pickle.dump(database95, file)
file.close()  

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
%reset -f
%matplotlib qt

import os
import pickle
import numpy as np
import itertools
import time
import pandas as pd
import math

#%% OPEN 96_ANTIBIOTICS

# A1/2/3 + B0 
os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\96_Antibiotics_A123+B0\\"

file = open(path + "/db_m_part_pro.pkl", "rb")
db_m_part_pro = pickle.load(file)

file = open(path + "/dict_model_pro.pkl", "rb")
dict_model_pro = pickle.load(file)

#%%

"""
96_Antibiotics contains all data of A123+B0 with antibiotic concentrations of 5, 5.75, 7.0 and 20. 
96_Antibiotics is added to database93, which already consits the data of A123+B0 with antibiotic concentration of 0.
Note: not all virtype combinations are simulated with antibiotics, so some places will be empty.
"""

#%%

os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\00_Complete_Database\\"

file = open(path + "/database93.pkl", "rb")
database93 = pickle.load(file)

#%% ADD MEANS TO DATABASE
# from dict_model_pro into database 'mean'

for virtype in dict_model_pro:
    for virstrength_abconc in dict_model_pro[virtype]:
        virstrength = virstrength_abconc[:-1]
        abconc = virstrength_abconc[-1]
        
        database93[virtype][virstrength][abconc]['mean'] = dict_model_pro[virtype][virstrength_abconc]
        
#%% ADD ITERATIONS TO DATABASE
# from dict_model_pro into database 'iterations'

for virtype in db_m_part_pro:
    for virstrength_abconc_it in db_m_part_pro[virtype]:
        virstrength = virstrength_abconc_it[:2]
        abconc = virstrength_abconc_it[2]
        it = virstrength_abconc_it[-1]
        
        database93[virtype][virstrength][abconc]['iterations'][it] = db_m_part_pro[virtype][virstrength_abconc_it]

#%% SAFE FILE - DATABASE 93

os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\00_Complete_Database\\"

file = open(path + "/database93.pkl", "wb")
pickle.dump(database93, file)
file.close()  

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
%reset -f
%matplotlib qt

import os
import pickle
import numpy as np
import itertools
import time
import pandas as pd
import math

#%% OPEN 97_ANTIBIOTICS

# A1/2/3 + B0 
os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\97_Antibiotics_A1+B1\\"

file = open(path + "/db_m_part_pro.pkl", "rb")
db_m_part_pro = pickle.load(file)

file = open(path + "/dict_model_pro.pkl", "rb")
dict_model_pro = pickle.load(file)

#%%

"""
97_Antibiotics contains all data of A1+B1 with antibiotic concentrations of 5, 5.75, 7.0 and 20. 
97_Antibiotics is added to database94, which already consits the data of A1+B1 with antibiotic concentration of 0.
Note: not all virtype combinations are simulated with antibiotics, so some places will be empty.
"""

#%%

os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\00_Complete_Database\\"

file = open(path + "/database94.pkl", "rb")
database94 = pickle.load(file)

#%%

"""
In 97_Antibiotics, combinations of vC and vS was switched. So instead of Ac + Bs, the simulation was run with As + Bc. 
Therefor this first needs to be switched again before this data can be added to database94. 
""" 

#%% SWITCH 'As + Bc' TO 'Ac + Bs' IN DICT_MODEL_PRO

# create new key for Ac + Bs
dict_model_pro['Ac + Bs'] = {}

for key in dict_model_pro['As + Bc']:
    new_key = (key[1],) + (key[0],) + (key[-1],)
    
    dict_model_pro['Ac + Bs'][new_key] = {}
    dict_model_pro['Ac + Bs'][new_key] = dict_model_pro['As + Bc'][key].copy()

new_column_names = {'m-#A': 'm-#B',
                    'm-#B': 'm-#A',
                    'm-slope-#A': 'm-slope-#B',
                    'm-slope-#B': 'm-slope-#A',
                    'm-m-Phag_prob_M_A': 'm-m-Phag_prob_M_B', 
                    'm-m-Phag_prob_M_B': 'm-m-Phag_prob_M_A', 
                    'std-#A': 'std-#B', 
                    'std-#B': 'std-#A', 
                    'std-slope-#A': 'std-slope-#B', 
                    'std-slope-#B': 'std-slope-#A',
                    'std-m-Phag_prob_M_A': 'std-m-Phag_prob_M_B', 
                    'std-m-Phag_prob_M_B': 'std-m-Phag_prob_M_A'
                    }

column_order = ['Time', 'm-Number_microbes', 'm-#A', 'm-#B', 'm-slope-#A', 'm-slope-#B',
                'm-m-Phag_prob_M_A', 'm-m-Phag_prob_M_B', 'm-Number_neutrophils',
                'm-Number_recruited', 'm-m-Phagocytose_rate', 'm-Bacterial_survival',
                'std-Number_microbes', 'std-#A', 'std-#B', 'std-slope-#A', 'std-slope-#B',
                'std-m-Phag_prob_M_A', 'std-m-Phag_prob_M_B', 'std-Number_neutrophils',
                'std-Number_recruited', 'std-m-Phagocytose_rate', 'std-Bacterial_survival', 'AB_fix_I']

for new_key in dict_model_pro['Ac + Bs']:
    # switch A and B
    dict_model_pro['Ac + Bs'][new_key].rename(columns=new_column_names, inplace=True)
    
    # re-order columns
    dict_model_pro['Ac + Bs'][new_key] = dict_model_pro['Ac + Bs'][new_key].reindex(columns=column_order)

# delete As + Bc    
del dict_model_pro['As + Bc']

#%% SWITCH 'As + Bc' TO 'Ac + Bs' IN DB_M_PART

# create new key for Ac + Bs
db_m_part_pro['Ac + Bs'] = {}

for key in db_m_part_pro['As + Bc']:
    new_key = (key[1],) + (key[0],) + (key[2],) + (key[3],)
    
    db_m_part_pro['Ac + Bs'][new_key] = {}
    db_m_part_pro['Ac + Bs'][new_key] = db_m_part_pro['As + Bc'][key].copy()

new_column_names = {'#A': '#B', 
                    '#B': '#A', 
                    'slope-#A': 'slope-#B', 
                    'slope-#B': 'slope-#A',
                    'm-Growthrate_A': 'm-Growthrate_B', 
                    'std-Growthrate_A': 'std-Growthrate_B',
                    'm-Growthrate_B': 'm-Growthrate_A', 
                    'std-Growthrate_B': 'std-Growthrate_A', 
                    'm-Killrate_A': 'm-Killrate_B', 
                    'std-Killrate_A': 'std-Killrate_B',
                    'm-Killrate_B': 'm-Killrate_A', 
                    'std-Killrate_B': 'std-Killrate_A', 
                    'm-NetGrowthrate_A': 'm-NetGrowthrate_B', 
                    'std-NetGrowthrate_A': 'std-NetGrowthrate_B',
                    'm-NetGrowthrate_B': 'm-NetGrowthrate_A', 
                    'std-NetGrowthrate_B': 'std-NetGrowthrate_A', 
                    'm-Phag_prob_M_A': 'm-Phag_prob_M_B',
                    'std-Phag_prob_M_A': 'std-Phag_prob_M_B', 
                    'm-Phag_prob_M_B': 'm-Phag_prob_M_A', 
                    'std-Phag_prob_M_B': 'std-Phag_prob_M_A' 
                    }

column_order = ['timestep', 'Time', 'Number_microbes', 'Number_neutrophils',
                'Number_recruited', 'AB_fix_I', '#A', '#B', 'slope-#A', 'slope-#B',
                'Bacterial_survival', 'iteration', 'm-Growthrate_A', 'std-Growthrate_A',
                'm-Growthrate_B', 'std-Growthrate_B', 'm-Killrate_A', 'std-Killrate_A',
                'm-Killrate_B', 'std-Killrate_B', 'm-NetGrowthrate_A', 'std-NetGrowthrate_A',
                'm-NetGrowthrate_B', 'std-NetGrowthrate_B', 'm-Phag_prob_M_A',
                'std-Phag_prob_M_A', 'm-Phag_prob_M_B', 'std-Phag_prob_M_B',
                'm-Number_phagocytosed', 'std-Number_phagocytosed', 'm-Phag_prob_N',
                'std-Phag_prob_N', 'm-Phagocytose_rate', 'std-Phagocytose_rate']

for new_key in db_m_part_pro['Ac + Bs']:
    # switch A and B
    db_m_part_pro['Ac + Bs'][new_key].rename(columns=new_column_names, inplace=True)
    
    # re-order columns
    db_m_part_pro['Ac + Bs'][new_key] = db_m_part_pro['Ac + Bs'][new_key].reindex(columns=column_order)

# delete As + Bc    
del db_m_part_pro['As + Bc']


#%% ADD MEANS TO DATABASE
# from dict_model_pro into database 'mean'

for virtype in dict_model_pro:
    for virstrength_abconc in dict_model_pro[virtype]:
        virstrength = virstrength_abconc[:-1]
        abconc = virstrength_abconc[-1]
        
        database94[virtype][virstrength][abconc]['mean'] = dict_model_pro[virtype][virstrength_abconc]
        
#%% ADD ITERATIONS TO DATABASE
# from dict_model_pro into database 'iterations'

for virtype in db_m_part_pro:
    for virstrength_abconc_it in db_m_part_pro[virtype]:
        virstrength = virstrength_abconc_it[:2]
        abconc = virstrength_abconc_it[2]
        it = virstrength_abconc_it[-1]
        
        database94[virtype][virstrength][abconc]['iterations'][it] = db_m_part_pro[virtype][virstrength_abconc_it]

#%% SAFE FILE - DATABASE 94

os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\00_Complete_Database\\"

file = open(path + "/database94.pkl", "wb")
pickle.dump(database94, file)
file.close()  

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
%reset -f
%matplotlib qt

import os
import pickle
import numpy as np
import itertools
import time
import pandas as pd
import math

#%% OPEN 98_ANTIBIOTICS

# A2 + B1
os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\98_Antibiotics_A2+B1\\"

file = open(path + "/db_m_part_pro.pkl", "rb")
db_m_part_pro = pickle.load(file)

file = open(path + "/dict_model_pro.pkl", "rb")
dict_model_pro = pickle.load(file)

#%%

"""
98_Antibiotics contains all data of A2+B1 with antibiotic concentrations of 5, 5.75, 7.0 and 20. 
98_Antibiotics is added to database95, which already consits the data of A2+B1 with antibiotic concentration of 0.
Note: not all virtype combinations are simulated with antibiotics, so some places will be empty.
"""

#%%

os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\00_Complete_Database\\"

file = open(path + "/database95.pkl", "rb")
database95 = pickle.load(file)

#%% ADD MEANS TO DATABASE
# from dict_model_pro into database 'mean'

for virtype in dict_model_pro:
    for virstrength_abconc in dict_model_pro[virtype]:
        virstrength = virstrength_abconc[:-1]
        abconc = virstrength_abconc[-1]
        
        database95[virtype][virstrength][abconc]['mean'] = dict_model_pro[virtype][virstrength_abconc]
        
#%% ADD ITERATIONS TO DATABASE
# from dict_model_pro into database 'iterations'

for virtype in db_m_part_pro:
    for virstrength_abconc_it in db_m_part_pro[virtype]:
        virstrength = virstrength_abconc_it[:2]
        abconc = virstrength_abconc_it[2]
        it = virstrength_abconc_it[-1]
        
        database95[virtype][virstrength][abconc]['iterations'][it] = db_m_part_pro[virtype][virstrength_abconc_it]

#%% SAFE FILE - DATABASE 95

os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\00_Complete_Database\\"

file = open(path + "/database95.pkl", "wb")
pickle.dump(database95, file)
file.close()  

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
%reset -f
%matplotlib qt

import os
import pickle
import numpy as np
import itertools
import time
import pandas as pd
import math

#%% OPEN 98_ANTIBIOTICS

os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\00_Complete_Database\\"

file = open(path + "/database93.pkl", "rb")
database93 = pickle.load(file)

file = open(path + "/database94.pkl", "rb")
database94 = pickle.load(file)

file = open(path + "/database95.pkl", "rb")
database95 = pickle.load(file)

#%% MERDE DATABASE93, DATABASE94 AND DATABASE95 INTO COMPLETE DATABASE

database = database93 | database94 | database95

#%% SAVE COMPLETE DATABASE

os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\00_Complete_Database\\"

file = open(path + "/database.pkl", "wb")
pickle.dump(database, file)
file.close() 

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

"""
Add settings 
"""

#%%
%reset -f
%matplotlib qt

import os
import pickle
import numpy as np
import itertools
import time
import pandas as pd
import math

#%% OPEN COMPLETE DATABASE

os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\00_Complete_Database\\"

file = open(path + "/database.pkl", "rb")
database = pickle.load(file)

#%% DEF FUNCTIONS
        
def open_all_settings(folder_path, folder_name):
    
    # change to correct folder
    os.chdir(r'D:\Research Project 1')
    folder_name = folder_name  
    path = os.getcwd() + "\\" + folder_path + "\\" + folder_name
    
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

abconcs = [5, 5.75, 7, 20]

# [[folder_path, folder_name, all_settings, virtype, abconc, virstrength]]

datas = [["03_Virulence_Strengths_A1+B0", "220612_vf_cellbound", 'Ac + B_', 0, None],
         ["03_Virulence_Strengths_A1+B0", "220612_vf_secreted", 'As + B_', 0, None],
         ["03_Virulence_Strengths_A1+B0", "220612_vf_recruit", 'Ar + B_', 0, None],
         
         ["04_Virulence_A2+B0", "220623_Acr_B", 'Acr + B_', 0, None],
         ["04_Virulence_A2+B0", "220623_Acs_B", 'Acs + B_', 0, None],
         ["04_Virulence_A2+B0", "220623_Asr_B", 'Asr + B_', 0, None],
         
         ["05_Virulence_A3+B0", "220623_Acsr_B", 'Acsr + B_', 0, None],
         
         ["06_Virulence_A1+B1", "220623_Ac_Br", 'Ac + Br', 0, None],
         ["06_Virulence_A1+B1", "220623_Ac_Bs", 'Ac + Bs', 0, None],
         ["06_Virulence_A1+B1", "220623_As_Br", 'As + Br', 0, None],
         
         ["07_Virulence_A2+B1", "220701_Acr_Bc", 'Acr + Bc', 0, None],
         ["07_Virulence_A2+B1", "220701_Acr_Br", 'Acr + Br', 0, None],
         ["07_Virulence_A2+B1", "220701_Acr_Bs", 'Acr + Bs', 0, None],
         
         ["07_Virulence_A2+B1", "220701_Acs_Bc", 'Acs + Bc', 0, None],
         ["07_Virulence_A2+B1", "220701_Acs_Br", 'Acs + Br', 0, None],
         ["07_Virulence_A2+B1", "220701_Acs_Bs", 'Acs + Bs', 0, None],
         
         ["07_Virulence_A2+B1", "220701_Ars_Bc", 'Asr + Bc', 0, None],
         ["07_Virulence_A2+B1", "220701_Ars_Br", 'Asr + Br', 0, None],
         ["07_Virulence_A2+B1", "220701_Ars_Bs", 'Asr + Bs', 0, None],
         
         ["20_Antibiotic_A123+B0", "", 'As + B_', abconcs, ('High', 'None')],
         ["20_Antibiotic_A123+B0", "", 'Asr + B_', abconcs, ('High', 'None')],
         ["20_Antibiotic_A123+B0", "", 'Acs + B_', abconcs, None],
         ["20_Antibiotic_A123+B0", "", 'Acsr + B_', abconcs, None],
         
         ["21_Antibiotic_A1+B1", "", 'Ac + Bs', abconcs, ('High', 'Low')],
         ["21_Antibiotic_A1+B1", "", 'Ac + Bs', abconcs, ('High', 'Medium')],
         ["21_Antibiotic_A1+B1", "", 'Ac + Bs', abconcs, ('High', 'High')],
         
         ["21_Antibiotic_A1+B1", "", 'As + Br', abconcs, ('Low', 'High')],
         ["21_Antibiotic_A1+B1", "", 'As + Br', abconcs, ('Medium', 'High')],
         ["21_Antibiotic_A1+B1", "", 'As + Br', abconcs, ('High', 'Medium')],
         ["21_Antibiotic_A1+B1", "", 'As + Br', abconcs, ('High', 'High')],
         
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Bc\\part1", 'Acs + Bc', abconcs, None],
         
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Bc\\part2", 'Acr + Bc', abconcs, ('Medium', 'Medium')],
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Bc\\part2", 'Acr + Bc', abconcs, ('High', 'High')],
         
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Bc\\part2", 'Asr + Bc', abconcs, ('Low', 'High')],
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Bc\\part2", 'Asr + Bc', abconcs, ('Medium', 'Medium')],
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Bc\\part2", 'Asr + Bc', abconcs, ('Medium', 'High')],
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Bc\\part2", 'Asr + Bc', abconcs, ('High', 'High')],
         
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Br\\part1", 'Acs + Br', abconcs, ('Low', 'Medium')],
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Br\\part1", 'Acs + Br', abconcs, ('Medium', 'Medium')],
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Br\\part1", 'Acs + Br', abconcs, ('High', 'Medium')],
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Br\\part1", 'Acs + Br', abconcs, ('Low', 'High')],
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Br\\part1", 'Acs + Br', abconcs, ('Medium', 'High')],
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Br\\part1", 'Acs + Br', abconcs, ('High', 'High')],
         
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Br\\part2", 'Asr + Br', abconcs, ('Low', 'High')],
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Br\\part2", 'Asr + Br', abconcs, ('Medium', 'High')],
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Br\\part2", 'Asr + Br', abconcs, ('High', 'Medium')],
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Br\\part2", 'Asr + Br', abconcs, ('High', 'High')],
         
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Bs\\part1", 'Acs + Bs', abconcs, None],
         
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Bs\\part2", 'Acr + Bs', abconcs, ('High', 'Low')],
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Bs\\part2", 'Acr + Bs', abconcs, ('High', 'Medium')],
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Bs\\part2", 'Acr + Bs', abconcs, ('High', 'High')],
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Bs\\part2", 'Acr + Bs', abconcs, ('Medium', 'Medium')],
         
         ["22_Antibiotic_A2+B1", "Antibiotics_A2_Bs\\part3", 'Asr + Bs', abconcs, None],
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
                
#%% clean up database, i.e. when no simulation is run, delete the empty dicts

del_keys = []
for virtype in database:
    for virstrength in database[virtype]:
        for abconc in database[virtype][virstrength]:            
            if len(database[virtype][virstrength][abconc]['mean']) == 0: # if empty  
                booleans = []
                if len(database[virtype][virstrength][abconc]['mean']) == 0:
                    booleans.append(False)
                if len(database[virtype][virstrength][abconc]['iterations']) == 0:
                    booleans.append(False)
                if len(database[virtype][virstrength][abconc]['settings']) == 0:
                    booleans.append(False)
                
                if len(booleans) != 3:
                    raise Exception('ERROR')
                
                del_keys.append([virtype, virstrength, abconc])
                
for key in del_keys:
    del database[key[0]][key[1]][key[2]]

#%% SAVE COMPLETE DATABASE INCL SETTINGS

os.chdir(r'D:\Research Project 1')
path = os.getcwd() + "\\00_Complete_Database\\"

file = open(path + "/database.pkl", "wb")
pickle.dump(database, file)
file.close() 

#%%

os.chdir(r'C:/Users/eske/Documents/Bio-Pharmaceutical Sciences/Jaar 1/Research Project 1/Database')
path = os.getcwd() + "\\00_Complete_Database\\"

file = open(path + "/database.pkl", "wb")
pickle.dump(database, file)
file.close() 













    

          
           
