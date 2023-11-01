# -*- coding: utf-8 -*-
"""
An agent based model of a polymicrobial community with virulence-mediated interspecies interactions
Author: Eske van Meegen and Cathi Herzberg, LACDR, Leiden University, The Netherlands

Parts of the code are based on Herzberg C. (2023) Vanhasseltlab/polymicro-abm: Polymicro-abm. v1.0.0. An agent
based modeling framework of an interactive, polymicrobial community during antimicrobial treatment. Zenodo. 
doi:10.5281/zenodo.10039816 which can be found on https://github.com/vanhasseltlab/PolyMicro-ABM 

File Description: PostSimulation script for Eske_Simulations_Model_v5 and Eske_Simulations_Model_v6. This script processes raw simulation data immediately after run.
"""
#%%
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------                              
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
##%% RELOAD RAW SIMULATION DATA
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------                             
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 

import pickle
import numpy as np
import pandas as pd
import time
import warnings

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------                             
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 
##%% POST-SIMULATION PROCESSING OF RAW DATA
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------                              
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- 


def main(data_interval,path,db_agents,db_model,all_species,ab_conc_range,microbe_capacity,n_init,neutrophil_capacity,neutrophil_n_init,its,runtime,gridrange,stepsize,data_collection,variable_params,fixed_params):

# ------------------------------------------------------------------------------------------------------------------     
# SPLIT DB_AGENTS INTO DB_MICROBES AND DB_NEUTROPHILS
# ------------------------------------------------------------------------------------------------------------------ 
    
    t0=time.time() # start timer   

    db_microbes = {}
    db_neutrophils = {}
    for key in db_agents:
        db_microbes[key] = db_agents[key][["Species", "Growthrate", "Killrate","NetGrowthrate","Phag_prob_M"]].dropna(how='all')       
        db_neutrophils[key] = db_agents[key][["Recruit_time","Age","Death_time","Number_phagocytosed","Phag_prob_N"]].dropna(how='all')
    
    db_agents = db_microbes # for now to not have to change every db_agents in the script to db_microbes from now on
    
    file = open(path + "/db_neutrophils.pkl", "wb")
    pickle.dump(db_neutrophils, file)
    file.close()  
    
    t1=time.time()
    #print("SPLIT DB_AGENTS: ",t1-t0)

# ------------------------------------------------------------------------------------------------------------------ 
# PROCESSING OF DB_AGENTS, DB_NEUTROPHILS AND DB_MODEL
# ------------------------------------------------------------------------------------------------------------------                              
    
    t0=time.time() # start timer
    
    
    # PROCES DB_MODEL
    
    # add columns, do calculations before reformatting
    for key in db_agents:       
        # split AB_fix column, multiple drugs
        for y in range(len(ab_conc_range[0])): # for every drug in the simulation
            db_model[key] = pd.concat([db_model[key],
                                       pd.DataFrame({'AB_fix_'+"I"*(y+1): [x[y] for x in db_model[key]["AB_fix"]]},
                                                    index=db_model[key].index)],
                                      axis=1, join='inner')
                     
        # remove old column
        db_model[key]=db_model[key].drop(['AB_fix'],axis=1)
    
    species_names=[x[0] for x in all_species]
    for key in db_model:
        
        # process number of bacteria per species
        for y in range(len(species_names)):
            db_model[key] = pd.concat([db_model[key],
                                       pd.DataFrame({'#'+species_names[y]: [x[species_names[y]] for x in db_model[key]['Number_species']]},index=db_model[key].index)],
                                      axis=1, join='inner')
        db_model[key]=db_model[key].drop(['Number_species'],axis=1) #remove old column
        
        # add column for slope of species growth
        for y in range(len(species_names)):
            neto = pd.DataFrame({"NET": db_model[key]["#"+species_names[y]].diff()})
            time_diff = pd.DataFrame({"time_diff": db_model[key]["Time"].diff()})
            
            db_model[key] = pd.concat([db_model[key],
                                       pd.DataFrame({"slope-#"+species_names[y]: neto["NET"]/time_diff["time_diff"]},
                                                    index=db_model[key].index)],
                                      axis=1, join='inner')        
        
        # add bacterial survival column for every iteration         
        bacterial_survival = [100] * len(db_model[key]) # 100 represents survival, i.e. there are still bacteria at the end of the simulation   
        number_microbes = list(db_model[key]["Number_microbes"])
        if number_microbes[-1] == 0: # if no bacteria are alive at the end of the simulation
            bacterial_survival[-1] = 0
      
        db_model[key] = pd.concat([db_model[key],
                                   pd.DataFrame({'Bacterial_survival': bacterial_survival},
                                                index = db_model[key].index)],
                                  axis = 1, join = 'inner') 
        
    # PROCES DB_AGENTS
    
    # add timestep, agentID and Time column to each dataframe in agent dictionary
    for key in db_agents:
        timestep = [x[0] for x in db_agents[key].index]
        time_list = list(db_model[key]['Time'])
        # generate Time column for agent dataframes
        time_a = []
        j=0
        timestep_0=timestep[0]
        for i in range(len(timestep)):
            if timestep[i]!=timestep_0:
                timestep_0=timestep[i]
                j+=1
            time_a.append(time_list[j])
        # add time and timepoint columns
        db_agents[key] = pd.concat([db_agents[key],
                                    pd.DataFrame({'timestep': timestep,
                                                  'AgentID':[x[1] for x in db_agents[key].index],
                                                  'Time': time_a},
                                                 index=db_agents[key].index)],
                                   axis=1, join='inner')
  
    # PROCES DB_NEUTROPHILS
    
    # add timestep, agentID and Time column to each dataframe in neutrophil dictionary
    for key in db_neutrophils:
        timestep = [x[0] for x in db_neutrophils[key].index] 
        time_list = list(db_model[key]['Time'])
        # generate Time column for neutrophil dataframes 
        time_n = []
        j = 0
        timestep_0 = timestep[0]
        for i in range(len(timestep)):
            if timestep[i] != timestep_0:
                timestep_0 = timestep[i]
                j += 1
            time_n.append(time_list[j]) # time_n is list of time at each step for each neutrophil agent
        # add timestep and time columns
        db_neutrophils[key] = pd.concat([db_neutrophils[key],
                                         pd.DataFrame({'timestep': timestep, 
                                                       'AgentID': [x[1] for x in db_neutrophils[key].index],
                                                       'Time': time_n},
                                                      index = db_neutrophils[key].index)], 
                                        axis=1, join='inner')
    
    # calculate number phagocytosed per min and add to each dataframe in neutrophil dictionary
    for key in db_neutrophils:  
        phag_diff = pd.DataFrame({"phag_diff": db_neutrophils[key].groupby(axis=0, level="AgentID")["Number_phagocytosed"].diff()})
        time_diff = pd.DataFrame({"time_diff": db_neutrophils[key].groupby(axis=0, level="AgentID")["Time"].diff()})
        
        for i,row in time_diff.iterrows():
            # if neutrophil is recruited inbetween data collection intervals, time and phag diff have to be calculated manually
            if np.isnan(time_diff["time_diff"].loc[i]) and not db_neutrophils[key]["Recruit_time"].loc[i] % data_interval == 0: 
                time_diff["time_diff"].loc[i] = db_neutrophils[key]["Time"].loc[i] - db_neutrophils[key]["Recruit_time"].loc[i]
                phag_diff["phag_diff"].loc[i] = db_neutrophils[key]["Number_phagocytosed"].loc[i]            
            # if neutrophil died inbetween data collection intervals, time_diff has to be calculated manually
            if not np.isnan(db_neutrophils[key]["Death_time"].loc[i]) and not db_neutrophils[key]["Death_time"].loc[i] % data_interval == 0: 
                time_diff.loc[i] = db_neutrophils[key]["Death_time"].loc[i] % data_interval
    
        # concatenate dataframe of phagocytose rate (=phag_diff/time_diff) to db_neutrophils
        db_neutrophils[key] = pd.concat([db_neutrophils[key],
                                         pd.DataFrame({"Phagocytose_rate": phag_diff["phag_diff"]/time_diff["time_diff"]})],
                                        axis=1, join='inner')
        
        
    # save processed version of raw data, all timepoints               
    file = open(path + "/db_agents-pro.pkl", "wb")
    pickle.dump(db_agents, file)
    file.close() 
    
    file = open(path + "/db_neutrophils-pro.pkl", "wb")
    pickle.dump(db_neutrophils, file)
    file.close()
    
    file = open(path + "/db_model-pro.pkl", "wb")
    pickle.dump(db_model, file)
    file.close()    
    
    t1=time.time()
    #print("DB_AGENTS, DB_NEUTROPHILS, DB_MODEL: ",t1-t0)

# ------------------------------------------------------------------------------------------------------------------ 
# DB_A_PART, DB_N_PART AND DB_M_PART
# ------------------------------------------------------------------------------------------------------------------ 
    
    t0 = time.time() # start timer

    # select part of the datafarmes for calculations and count number of each species
                   
    #DEFINE INTERVAL TO COLLECT SUBSET OF DATA FOR CALCULATIONS
    interval = data_interval # calculate species count etc every x mins, choose so that this interval >= datacollection interval
    t_endpoints_time =  list(range(0,round(runtime),interval)) + [runtime]
    
    # determine species names for all simulated scenarios
    species_names=[x[0] for x in all_species]
    
    # take a subset of the dataframe at the defined endpoints
    db_m_part = db_model.copy()
    db_a_part = db_agents.copy()
    db_n_part = db_neutrophils.copy()
    del [db_model, db_agents, db_neutrophils]
   
    for key in db_m_part:
        time_list = list(db_m_part[key]["Time"])
        Number_microbes=list(db_m_part[key]["Number_microbes"])
        
        t_endpoints_a = [x for x in t_endpoints_time if x<=time_list[-1]+interval] # remove endpoints that haven't been simulated
        mins = [[abs(x-y) for y in time_list] for x in t_endpoints_a]
        min_index = [i for x in mins for i in range(len(x)) if min(x)==x[i]] # index of which timestep is closest to endpoint
        t_endpoints_a = [time_list[i] for i in min_index] # generate endpoints from timesteps closest
    
        db_m_part[key] = db_m_part[key][db_m_part[key]['Time'].isin(t_endpoints_a)]
        db_a_part[key] = db_a_part[key][db_a_part[key]['Time'].isin(t_endpoints_a)]
        db_n_part[key] = db_n_part[key][db_n_part[key]['Time'].isin(t_endpoints_a)]

    # add iteration number to row
    for key in db_m_part:
        l = len(db_m_part[key])
        db_m_part[key]=pd.concat([db_m_part[key],
                                  pd.DataFrame({'iteration' :[key[-1]]*l},
                                               index=db_m_part[key].index)], 
                                 axis=1)
    
    # caclulate mean and std for some columns from agent dictionary
    # this is only interesting if variable in column are heterogeneous across population
    column_names = ['Growthrate','Killrate',"NetGrowthrate","Phag_prob_M"]#+['MIC_'+'I'*(y+1) for y in range(len(ab_conc_range[0]))]
         
    for key in db_m_part:
        tp=list(dict.fromkeys([x[0] for x in db_a_part[key].index])) # timepoints in agent dataframe
        indexes=[]
        for t in tp: # go through all time points
            # make lists with indexes per species for each time point
            indexes_species = [[tuple([t,index]) for index,row in db_a_part[key].loc[t].iterrows() if row["Species"]==spec_name] for spec_name in species_names]
            indexes.append(indexes_species)      
        indexes = [[x[i] for x in indexes] for i in range(len(indexes[0]))] # resort, switch 1st and 2nd axis
            
        for col in column_names:
            for s in range(len(species_names)):
                # calculate mean and std for each column in list
                mn = [np.mean([db_a_part[key].loc[x_t][col]]) for x_t in indexes[s]]
                stds = [np.std([db_a_part[key].loc[x_t][col]]) for x_t in indexes[s]]
                
                # if last timepoint doesn't exist in db_a_part, because all microbes are dead
                Number_microbes = list(db_m_part[key]["Number_microbes"])
                if Number_microbes[-1] == 0:
                    mn.append(np.nan)
                    stds.append(np.nan)
                
                # add to db_m_part
                db_m_part[key]=pd.concat([db_m_part[key],
                                          pd.DataFrame({'m-'+col + '_'+ str(species_names[s]) :mn ,
                                                        'std-'+col + '_'+ str(species_names[s]) :stds},
                                                       index=db_m_part[key].index)],
                                         axis=1, join='inner')               
                
                
    # calculate mean and std for some columns from neutrophil dictionary
    # add mean and std column to db_m_part for selected columns per timepoint
    np_column_names = ['Number_phagocytosed','Phag_prob_N','Phagocytose_rate']  

    for key in db_m_part:
        tp = list(dict.fromkeys([x[0] for x in db_n_part[key].index])) # timepoints in neutrophil dataframe
        all_steps = [] # list of dataframes of columns per timepoint per agent
        for t in tp:
            step = db_n_part[key].loc[t][np_column_names] # columns per timepoint per agent
            all_steps.append(step)
        
        for col in np_column_names:            
            mn = [np.mean(step[col]) for step in all_steps]
            std = [np.std(step[col]) for step in all_steps]
            
            db_m_part[key] = pd.merge(db_m_part[key], 
                                      pd.DataFrame({'m-'+col: mn,
                                                    'std-'+col: std}, 
                                                   index = tp), 
                                      left_on='timestep', right_index=True, how='left')
            
        
    # save processed data of only select timepoints
    file = open(path + "/db_m_part.pkl", "wb")
    pickle.dump(db_m_part, file)
    file.close()
    
    file = open(path + "/db_a_part.pkl", "wb")
    pickle.dump(db_a_part, file)
    file.close()    

    file = open(path + "/db_n_part.pkl", "wb")
    pickle.dump(db_n_part, file)
    file.close()
    
    t1=time.time()
    #print('DB_A_PART, DB_N_PART AND DB_M_PART: ',t1-t0)

# ------------------------------------------------------------------------------------------------------------------ 
# DB_PER_NEUTROPHIL
# ------------------------------------------------------------------------------------------------------------------  
# =============================================================================
#     t0 = time.time() # start timer    
#     
#     # creates a database which is a dict (keys: simulations) of dict (keys: agentIDs) of dataframes (data per neutrophil for all timepoints)
#     # in order to plot data of each seperate neutrophil agent for each simulation
#     # all timepoints of a simulation have to be in the dataframe of each agent otherwise it can't be plotted
# 
#     db_n_id = db_n_part.copy()
#     del [db_n_part]
# 
#     db_per_neutrophil = {} 
# 
#     for key in db_n_id: # each key represent a simulation
#         db_per_neutrophil[key] = {}
#     
#         agentIDs = list(set(db_n_id[key]["AgentID"])) # list of all agentIDs for this simulation
#         time_i = list(set(db_n_id[key]["Time"])) # list of all timepoints for this simulation
#         time_i.sort()
#     
#         # generate empty dataframe length of time_i
#         df_empty = pd.DataFrame({"Time": time_i})    
#     
#         # sort the dataframe of this simulation by the AgentID-index and remove the step-index
#         db_n_id[key] = db_n_id[key].sort_index(level=1).droplevel(0)
# 
#         for agentID in agentIDs:
#             df_neutro = db_n_id[key].loc[[agentID]].merge(df_empty, how='right', on="Time") # dataframe of data per neutrophil for all timepoints
#             db_per_neutrophil[key][agentID] = df_neutro 
#     
#     file = open(path + "/db_per_neutrophil.pkl", "wb")
#     pickle.dump(db_per_neutrophil, file)
#     file.close()
#     
#     t1=time.time()
#     #print('DB_PER_NEUTROPHIL: ',t1-t0)     
# =============================================================================
    
# ------------------------------------------------------------------------------------------------------------------ 
# DICT_MODEL
# ------------------------------------------------------------------------------------------------------------------ 

    t0 = time.time() # start timer1
    
    # Dictionary of Dataframe with info about all iterations together, average values etc. per combination: dict_model
    # Calculate mean and std of species number between iterations at each time point
 
    # make a copy only with select columns  
    column_s = ['iteration','Time','Number_microbes']+['#'+str(species_names[i]) for i in range(len(species_names))]+['slope-#'+str(species_names[i]) for i in range(len(species_names))]+['m-Phag_prob_M_'+str(species_names[i]) for i in range(len(species_names))]+['Number_neutrophils','Number_recruited','m-Phagocytose_rate','Bacterial_survival']+[]
    
    # make a dict (keys: simulations) of dataframes (indeces: steps) with selected columns
    db_m = db_m_part.copy()
    for key in db_m:
        db_m[key] = db_m[key][column_s]

    all_keys = list(db_m.keys()) # keys for all simulation (all abconc and its)
    keys = list(dict.fromkeys([x[:-1] for x in all_keys])) # key for each abconc (without its), the key list for the new dict_model

    dfs=[] # list of dataframes for the new dict_model

    for key in keys: # each key represents an abconc
        key_i=[key+(i,) for i in range(its)] # list of keys for each iteration of this abconc
        
        max_len = max([len(db_m[k_i]) for k_i in key_i]) # determine maximal length of simulation for this abconc
        last_time = [list(db_m[k_i]['Time'])[-1] for k_i in key_i] # make a list of last timepoints of simulations for this abconc
        same_lengths = all([len(db_m[k_i])==max_len for k_i in key_i]) # this is set to false when not all iteration have the same length
        same_times = all([list(db_m[k_i]['Time'])[-1]==last_time[0] for k_i in key_i]) # this is set to false when not all last time points are the same
        
        if not same_lengths:
            # make a df "empty" with the length of the longest iteration and only the index and Time column
            # concat to the shorter iterations a slice of this df from last timepoint of shorter iteration to the last timepoint of longest iteration
            key_longest_it = key_i[last_time.index(max(last_time))] # key of the longest iteration for this abconc
            df_empty = db_m[key_longest_it][[]] # dataframe with only index of longest iteration for this abconc
            
            for k_i in key_i:
                if len(db_m[k_i]) != max_len:
                    last_i = len(db_m[k_i])-1 # index of last timepoint
                    db_m[k_i] = pd.concat([db_m[k_i],df_empty[last_i+1:]]) # add df_empty from last index+1
                    db_m[k_i]['iteration'] = db_m[k_i]['iteration'].replace(np.nan, k_i[-1]) # add iteration number to newly added rows
                    db_m[k_i]['Bacterial_survival'] = db_m[k_i]['Bacterial_survival'].replace(np.nan, 0) # 0 represent no bacteria alive

# =============================================================================
#         if not same_lengths:
#             for k_i in key_i:
#                 while len(db_m[k_i]) is not max_len:
#                     db_m[k_i]=db_m[k_i].append(pd.DataFrame([[list(k_i)[-1]]+[np.nan,0]+[0]*len(species_names)+[np.nan]*len(species_names)+[np.nan,np.nan,np.nan,0]],columns=column_s),ignore_index=True) 
#                     # in the order of column_s
# =============================================================================
        
        same_lengths = all([len(db_m[k_i])==max_len for k_i in key_i]) # recalculate
        
        # for some of the col in column_s all values in the last row are np.nan
        # this raises a RuntimeWarning: mean of empty slice
        # because it is only the last row, it is not a concerning runtime issue
        with warnings.catch_warnings(): # this functions suppresses the warning during means calculation
            warnings.simplefilter("ignore", category=RuntimeWarning) 
        
            if same_lengths:
                mean_time_vector = [int(np.nanmean([x[i] for x in [list(db_m[k_i]['Time']) for k_i in key_i]])) for i in range(max_len)]
                means = [[np.nanmean([x[i] for x in [list(db_m[k_i][col]) for k_i in key_i]]) for i in range(max_len)] for col in column_s[2::]]
                stds = [[np.nanstd([x[i] for x in [list(db_m[k_i][col]) for k_i in key_i]]) for i in range(max_len)] for col in column_s[2::]]
                x=[mean_time_vector]+means+stds
                x=[[y[i] for y in x ] for i in range(len(x[0]))]
                columns_new=[column_s[1]]+['m-'+ x for x in column_s[2::]]+['std-'+ x for x in column_s[2::]]
            dm_means=pd.DataFrame(x,columns=columns_new)
            dfs.append(dm_means)

    del [db_m]
    dict_model = dict(zip(keys,dfs))
    
    # add AB_fix to dict_model
    for key in dict_model:
        for y in range(len(ab_conc_range[0])): # for every drug in the simulation
            ab_fix = db_m_part[key+(0,)]['AB_fix_'+'I'*(y+1)][0]
            df_ab_fix = pd.DataFrame({'AB_fix_'+'I'*(y+1): [ab_fix]*len(dict_model[key])},
                                     index = dict_model[key].index)
            dict_model[key] = pd.concat([dict_model[key],df_ab_fix],axis=1)  
    
           
    # save processed data with averages over all iterations
    file = open(path + "/dict_model.pkl", "wb")
    pickle.dump(dict_model, file)
    file.close()  
       
    t1=time.time()
    #print('DICT_MODEL: ',t1-t0)         
