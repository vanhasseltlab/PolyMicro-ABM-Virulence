# -*- coding: utf-8 -*-
"""
An agent based model of a polymicrobial community with virulence-mediated interspecies interactions
Author: Eske van Meegen, LACDR, Leiden University, The Netherlands

Parts of the code are based on Herzberg C. (2023) Vanhasseltlab/polymicro-abm: Polymicro-abm. v1.0.0. An agent
based modeling framework of an interactive, polymicrobial community during antimicrobial treatment. Zenodo. 
doi:10.5281/zenodo.10039816 which can be found on https://github.com/vanhasseltlab/PolyMicro-ABM 

File Description: Simulation script for a dualspecies population during antibiotic treatment.
"""

import os
import pickle
from mesa.batchrunner import BatchRunner
import numpy as np
import time
from Eske_Model_v6 import MyModel
from Eske_PostSim_v4 import main
#%% SIMULATION SETTINGS 

#-----------------------------------------------------------------------------------------------------------------------
# SET MICROBE SPECIES CHARACTERISTICS
#-----------------------------------------------------------------------------------------------------------------------

### SHARED CHARACTERISTICS OVER SPECIES ###
maxgrowthrate = 0.015 # 1/min , approximately 66 min generation time
species_MIC = [5] # list of MIC values of species A to every drug in the simulation, in this example there is only one drug, so it's a list with one element
k_d0 = 0.2/60 # natural death rate
phag_prob_M = 1 # probability microbe is going to get phagocytosed by a neutrophil upon interaction
phag_prob_M_min = 1/5000 #1/microbe_capacity # minimum phagocytose probability 

# characteristics of one species in the format [species name, species_MIC, maximum growth rate, maximum kill rate, phagocytose probabilities]
species_A = ["A", species_MIC, maxgrowthrate, k_d0, phag_prob_M, phag_prob_M_min] 
species_B = ["B", species_MIC, maxgrowthrate, k_d0, phag_prob_M, phag_prob_M_min]

all_species = [species_A, species_B] # all_species is a list of the characteristics of all species

#-----------------------------------------------------------------------------------------------------------------------
# SET MICROBE SPECIES VIRULENCE 
#-----------------------------------------------------------------------------------------------------------------------

virulence_types = ['cellbound', 'secreted', 'recruit']
virulence_strengths = ['Low', 'Medium', 'High']
virulence_strength_values = [{'None': 0, 'Low': 0.4, 'Medium': 0.8, 'High':0.9},           # cellbound
                             {'None': 0, 'Low': 0.025, 'Medium': 0.045, 'High':0.06},      # secreted
                             {'None': 0, 'Low': 0.5, 'Medium': 0.7, 'High':0.9}]           # recruit
 
all_species_all_virulence = [[['None', 'High', 'None'],['None', 'None', 'Low']]]

#%%  
#-----------------------------------------------------------------------------------------------------------------------
# SET NEUTROPHIL CHARACTERISTICS
#-----------------------------------------------------------------------------------------------------------------------

max_recruit_rate = 3/24/60 # 1/min #0.0045 
halflife = 7.19*60 # min
phag_rate = 0.081 # 1/min # phagocytose rate constant

neutrophils = [max_recruit_rate, halflife, phag_rate]

#-----------------------------------------------------------------------------------------------------------------------
# SET ANTIBIOTIC TREATMENT
#-----------------------------------------------------------------------------------------------------------------------

# drugs
k = [4] #[5] # hypothetical (original) & Tobramycin & Ceftazidime , hill parameters for each drug # Mouton& Vink Table 1 # independent of units
G_min = [None] #[-3/60] # for bactericidal antibiotics only
drugtype = ['S'] #['C'] # ADD static 'S' or for cidal 'C' to the drug list

# create list of all_drugs to be simulated in list 
all_drugs = [[[G_min[i]],[k[i]],[drugtype[i]]] for i in range(len(k))]
# these are four hypothetical antibiotics: bactericidal & Concentration-independent, bactericidal & Concentration-dependent, bacteriostatic & Concentration-independent, bacteriostatic & Concentration-dependent,

# concentrations simulated for each drug
#ab_conc_range_CI = [[x] for x in [x*species_MIC[0] for x in [0]]]#[1,1.3,1.65,4]]]
ab_conc_range_SI = [[x] for x in [x*species_MIC[0] for x in [0,1,1.15,1.4,4]]]
ab_conc_range_alldrugs = [ab_conc_range_SI]#,ab_conc_range_CI]

#-----------------------------------------------------------------------------------------------------------------------
# SET SIMULATION PARAMETERS
#-----------------------------------------------------------------------------------------------------------------------

# init and capacity
microbe_capacity = 5000 # maximal number of species in simulation at the same time
n_init = microbe_capacity # number of bacteria at the start of the simulation 
neutrophil_capacity = 100 # maximal number of species in simulation at the same time
neutrophil_n_init = 50 #int(0.5 * neutrophil_capacity) # number of neutrophils at the start of the simulation 
#n_n_inits = [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

# simulation settings
gridrange = 10 # width/length of square grid in nunber of patches, total number of patches = gridrange*gridrange
runtime = 1*24*60 # runtime in mins
its = 25 # number of iterations
data_collection = 20 # interval for datacollection in mins
data_agents = True # set to False when you don't want to collect agent reporters
stepsize = [5] # length of one timestep in simulation, in mins 
max_steps = round(runtime/stepsize[0]) # maximal number of steps to be performed

#print("Named simulation folder?")
#%% # BATCH RUN SIMULATIONS
# ------------------------------------------------------------------------------------------------------------------                              
# ------------------------------------------------------------------------------------------------------------------ 

time_start=time.time() # start timer

db_agents={}
db_model={}

# run simulations in batch 
# for one species at different drug concentrations

for drugs,ab_conc_range in zip(all_drugs,ab_conc_range_alldrugs):
    ab_conc_time_range = [[list(np.repeat(x, max_steps+1)) for x in ab_conc] for ab_conc in ab_conc_range] # create vectors with concentration over time at every timestep for each defined concentration
    for ab_conc,ab_conc_time in zip(ab_conc_range,ab_conc_time_range):
        for all_species_virulence in all_species_all_virulence:
          
            fixed_params = {
                            "microbe_capacity": microbe_capacity, "all_species": all_species, "all_species_virulence": all_species_virulence, "virulence_strength_values": virulence_strength_values,
                            "neutrophil_n": neutrophil_n_init, "neutrophil_capacity": neutrophil_capacity, "neutrophils": neutrophils,
                            "ab_conc_time": ab_conc_time, "drugs": drugs, 
                            "runtime": runtime, "gridrange": gridrange, "stepsize": stepsize[0], 
                            "data_collection": data_collection, "data_agents": data_agents
                            }
        
            variable_params = {"microbe_n": range(n_init,n_init+1)}

            batch_run = BatchRunner(MyModel,
                                    variable_params,
                                    fixed_params,
                                    iterations=its, # this defines that there are 3 iteratiins run for each scenario specified
                                    max_steps=max_steps)
            batch_run.run_all()
        
            #Get the Agent DataCollection
            data_batch_agents = batch_run.get_collector_agents()
            #Get the Model DataCollection
            data_batch_model = batch_run.get_collector_model()
            # every simulation in batch run will collect results in two dictionaries: data_batch_agents, data_batch_model
    
            # raw data from a batch run is stored in the dictionaries: db_agents, db_model
            # copy ditionaries to another dictionary with all the antibiotic concentrations, to do calculations across iterations
            j = -1
            for i in range(its):
                #print(drugs[1][0],ab_conc[0],i)
                j += 1
                if len(all_species) == 1:
                    db_agents[(ab_conc[0],
                               all_species_virulence[0][0],all_species_virulence[0][1],all_species_virulence[0][2],
                               i)]= data_batch_agents[(n_init,j)].copy()
                    
                    db_model[(ab_conc[0],
                              all_species_virulence[0][0],all_species_virulence[0][1],all_species_virulence[0][2],
                              i)]= data_batch_model[(n_init,j)].copy()
                    
                if len(all_species) == 2:
                    db_agents[(ab_conc[0],
                               all_species_virulence[0][0],all_species_virulence[0][1],all_species_virulence[0][2],
                               all_species_virulence[1][0],all_species_virulence[1][1],all_species_virulence[1][2],
                               i)]= data_batch_agents[(n_init,j)].copy()
                    
                    db_model[(ab_conc[0],
                              all_species_virulence[0][0],all_species_virulence[0][1],all_species_virulence[0][2],
                              all_species_virulence[1][0],all_species_virulence[1][1],all_species_virulence[1][2],
                              i)]= data_batch_model[(n_init,j)].copy()
    
            del [data_batch_agents,data_batch_model] # delete single dictionaries

# save diraw data to files  
os.chdir(r'C:\Users\eskev\Documents\Bio-Pharmaceutical Sciences - LU\Jaar 1\Research Project 1')
folder_name = "/21_Antibiotic_A1_B1"  
path = os.getcwd() + "\Simulations\\LastMin\\" + folder_name
os.mkdir(path)

file = open(path + "/db_agents.pkl", "wb")
pickle.dump(db_agents, file)
file.close()    

file = open(path + "/db_model.pkl", "wb")
pickle.dump(db_model, file)
file.close()   

file = open(path + "/species.pkl", "wb")
pickle.dump([all_species,all_species_all_virulence,ab_conc_range], file)
file.close()  

file = open(path + "/neutrophils.pkl", "wb")
pickle.dump([neutrophils], file)
file.close()  

file = open(path + "/settings.pkl", "wb")
pickle.dump([microbe_capacity,n_init,neutrophil_capacity,neutrophil_n_init,its,runtime,gridrange,stepsize,data_collection,variable_params,fixed_params], file)
file.close() 

time_end=time.time()
print("SIMULATION TOTAL TIME: ",str(time_end-time_start))


# ------------------------------------------------------------------------------------------------------------------                              
# ------------------------------------------------------------------------------------------------------------------ 

t0=time.time()

data_interval = data_collection
main(data_interval,path,db_agents,db_model,all_species,ab_conc_range,microbe_capacity,n_init,neutrophil_capacity,neutrophil_n_init,its,runtime,gridrange,stepsize,data_collection,variable_params,fixed_params)
# arguments for main function: data_interval for data processing, path to raw data, all variables saved in pickle files

t1=time.time()
print("POST_SIM TOTAL TIME: ",str(t1-t0))

#%% GET DBs IN VARIABLE EXPLORER

# =============================================================================
# folder_name="220708_test"  
# os.chdir(r'C:\Users\eskev\Documents\Bio-Pharmaceutical Sciences - LU\Jaar 1\Research Project 1')
# path = os.getcwd() + "\Simulations\\000000_Test\\" + folder_name
# 
# # =============================================================================
# # folder_name="220612_vf_cellbound"  
# # os.chdir(r'D:\Research Project 1')
# # path = os.getcwd() + "\\03_Virulence_Strengths_A1+B0\\" + folder_name
# # =============================================================================
# 
# #
# 
# #file = open(path + "/db_neutrophils-pro.pkl", "rb")
# #db_neutrophils = pickle.load(file)
# #file.close()  
# 
# #file = open(path + "/db_agents-pro.pkl", "rb")
# #db_agents = pickle.load(file)
# #file.close()
# 
# #file = open(path + "/db_model-pro.pkl", "rb")
# #db_model = pickle.load(file)   
# #file.close()
# 
# file = open(path + "/species.pkl", "rb")
# all_species,all_species_all_virulence,ab_conc_range = pickle.load(file)
# file.close()
# 
# file = open(path + "/neutrophils.pkl", "rb")
# neutrophils = pickle.load(file)
# file.close()
# 
# file = open(path + "/settings.pkl", "rb")
# microbe_capacity,n_init,neutrophil_capacity,neutrophil_n_init,its,runtime,gridrange,stepsize,data_collection,variable_params,fixed_params = pickle.load(file)
# file.close()
# 
# # processed data per iteration
# file = open(path + "/db_m_part.pkl", "rb")
# db_m_part = pickle.load(file)
# file.close()
# #
# file = open(path + "/db_a_part.pkl", "rb")
# db_a_part = pickle.load(file)
# file.close()
# 
# file = open(path + "/db_n_part.pkl", "rb")
# db_n_part = pickle.load(file)
# file.close()
# 
# # processed data per neutrophil per iteration
# #file = open(path + "/db_per_neutrophil.pkl", "rb")
# #db_per_neutrophil = pickle.load(file)
# #file.close()
#  
# # processed averaged data
# file = open(path + "/dict_model.pkl", "rb")
# dict_model = pickle.load(file)
# file.close()
# 
# =============================================================================
