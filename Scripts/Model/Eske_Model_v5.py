# -*- coding: utf-8 -*-
"""
An agent based model of a polymicrobial community with virulence-mediated interspecies interactions
Author: Eske van Meegen, LACDR, Leiden University, The Netherlands

Parts of the code are based on Herzberg C. (2023) Vanhasseltlab/polymicro-abm: Polymicro-abm. v1.0.0. An agent
based modeling framework of an interactive, polymicrobial community during antimicrobial treatment. Zenodo. 
doi:10.5281/zenodo.10039816 which can be found on https://github.com/vanhasseltlab/PolyMicro-ABM 

File Description: Main model file version 5 (with numeric input for virulence strengths)
"""
#%%
# load packages
from mesa import Agent, Model
from mesa.time import RandomActivation
from mesa.time import StagedActivation
from mesa.space import MultiGrid
from mesa.datacollection import DataCollector
import numpy as np
import copy

# agent object class
class MyMicrobe(Agent):
    
# AGENT INITIALIZATION     

    def __init__(self, model, microbe_id, species, virulence_factors):#, interactions, interactions_env):
        super().__init__(microbe_id, model)
        self.alive = True
        # Species specific parameters, these parameters do not change during the simulation
        self.species = species # variable containing most important species characteristics
        self.species_name = self.species[0] # name of species
        self.MIC_original = self.species[1] # original MIC for this species # turned into list for all drugs [MIC for drug I, MIC for drug II]
        self.k_repmax_original = self.species[2] # maximal specific growth rate for this species | original maximal replication rate for this species
        self.k_d0 = self.species[3]   

        # Parameters specific to each agent
        self.microbe_id = microbe_id # unique agent name
        self.MIC = copy.deepcopy(self.MIC_original) # current MIC of this agent, initialized to be species specific MIC # turned into list for all drugs, can be different from MIC_origina   
        self.k_repmax = self.k_repmax_original # current k_repmax for this agent, altered by interactions
        self.set_netgrowth(self.model.ab_conc) # this creates K_NET , replicationrate, deathrate, E_C and E_S for each time step

        # Species specific phagocytose and virulence parameters
        self.virulence_factors = virulence_factors
        self.vf_phag_surface = self.virulence_factors[0]
        self.vf_phag_excreted = self.virulence_factors[1]
        self.phag_prob_M_original = self.species[4] - self.vf_phag_surface # original baseline phagocytose probability
        self.phag_prob_M_min = self.species[5] # minimum phagocytose probability 

        # Phagocytose and virulence parameters/counters specific to each agent
        self.phag_prob_M = self.phag_prob_M_original # specific to each agent, can be altered by virulent microbes in patch
        self.number_vir_interactions = 0 # counter, resets every timestep
        
# STEP FUNCTION FOR MICROBES

    def step(self):
        self.move_agent()
        self.virulence_act()        
        self.die_divide()
    
# FUNCTIONS TO SET MICROBE CHARACTERISTICS
        
    def set_netgrowth(self, conc):
        k_rep = self.set_k_rep(self.k_repmax)
        self.set_drugeffect(conc)
        
        self.replicationrate = k_rep* self.E_S # actual replication rate including drug and interaction effect at this time point
        self.deathrate = self.k_d0* self.E_C # actual death rate including drug and interaction effect at this time point
        
        self.K_NET = self.replicationrate - self.deathrate # overall net growth rate including drug effect and interaction effects
           
    def set_k_rep(self,k_repmax):
        k_d0 = self.k_d0

        if self.model.capacity_type=='local' and self.model.replication_radius < (self.model.gridrange-1)/2: # if replication radius bigger than a certain value, it is actually global
            p = self.model.gridrange*self.model.gridrange # number of patches
            if self.model.replication_radius == 0:
                p_N = 1 #number of patches to be considered for replication capacity
                B = self.number_microbes_grid[self.pos[0]][self.pos[1]] # THIS IS THE VALUE OF b IF RADIUS =0 and only this patch is taken into account
            else:
                patches = self.model.grid.get_neighborhood(self.pos,moore=True,include_center=True,radius=self.model.replication_radius)
                p_N = len(patches)
                B = sum([self.model.number_microbes_grid[x[0]][x[1]] for x in patches])
            
            K =  self.model.capacity/p * p_N# capacity per patch [number of cells], multiplied by the number of patches in the capacity radius
        # ONLY USE GLOBAL SETTING
        elif self.model.capacity_type=='global' or self.model.replication_radius >= (self.model.gridrange-1)/2:
            B = self.model.number_microbes
            K = self.model.microbe_capacity
            
        k_rep = k_repmax *(1- B/K*(1- k_d0/k_repmax))
        return k_rep
    
    def set_drugeffect(self, conc):
        k_rep = self.set_k_rep(self.k_repmax_original) # replication rate disregarding interaction effects or drug effects
        all_E_S=[]
        all_E_C=[]
        
        for drugnumber in range(len(conc)):
            if self.model.drugtype[drugnumber] == 'C': # for bactericidal drugs
                k = self.model.k[drugnumber]
                G_min = self.model.G_min[drugnumber]
                MIC = self.MIC[drugnumber] # already including the interaction effects
                c = conc[drugnumber] # already including the interaction effects
                k_d0 = self.k_d0
                
                E_C = 1 + (((k_rep-G_min)/k_d0 - 1)*(c/MIC)**k) / ((c/MIC)**k + G_min/(k_d0-k_rep))
                E_S = 1 # no effect
                
            elif self.model.drugtype[drugnumber] == 'S':  # for bacteriostatic drugs
                k = self.model.k[drugnumber]
                MIC = self.MIC[drugnumber] # already including the interaction effects
                c = conc[drugnumber] # already including the interaction effects
                k_d0 = self.k_d0
                
                E_S = 1 - ((c/MIC)**k)/(((c/MIC)**k) + k_d0/(k_rep-k_d0))
                E_C = 1 # no effect
                
            all_E_S.append(E_S)
            all_E_C.append(E_C)
        
        # combine effects of multiple drugs at one time point using equation for drug combinations
        if sum([1 for x in conc if x>0])==0: # if all drug conc zero at this time point
            self.E_S = 1 # no drug effect
            self.E_C = 1 # no drug effect
        elif sum([1 for x in conc if x>0])==1: # if one drug is positive at this timepoint
            index = [i if all_E_S[i]!=1 else 0 for i in range(len(conc))][0] # E_S is between 0 and 1, 1 meaning there is no effect
            self.E_S = all_E_S[index]

            index = [i if all_E_C[i]!=1 else 0 for i in range(len(conc))][0] # E_C is between 1 and (k_rep-G_min)/k_d0, 1 meaning there is no effect
            self.E_C = all_E_C[index]

# ACTION FUNCTIONS FOR MICROBES

    # step function called at every time step for each agent by the scheduler
    def move_agent(self):
        if self.alive:
            for n in range(self.model.n_steps): # n_steps for each minute
#                self.model.number_microbes_grid[self.pos[0]][self.pos[1]] -=1 # remove counter from old position
                # 1) MOVE on grid
                if self.model.move_type == "stepwise":
                    self.move()
                elif self.model.move_type == "random":
                    self.move_rand()
#                self.model.number_microbes_grid[self.pos[0]][self.pos[1]] +=1 # add counter at new position            
     
    def move(self): # move to one of 8 neighbouring patches or stay on the same patch, all with equal chances
        possible_steps = self.model.grid.get_neighborhood(
            self.pos,
            moore=True,
            include_center=True)
        new_position = self.random.choice(possible_steps)
        self.model.grid.move_agent(self, new_position)    
    
    def virulence_act(self):
        self.phag_prob_M = self.phag_prob_M_original # reset phag_prob_M
        self.number_vir_interactions = 0 # reset counter for number of virulence interactions 
        
        # make list of all possible interaction agents
        microbe_cellmates = self.model.give_microbe_cellmates(self.pos) # this list also contains self, which is good because excreted molecule can also protect microbe that excreted the molecule
        for microbe in microbe_cellmates:
            if microbe.alive and microbe.vf_phag_excreted != 0: # if the microbe cellmate excretes virulence molecules 
                if self.phag_prob_M - microbe.vf_phag_excreted > 0: # phagocytose probability can't be smaller than 0
                    self.phag_prob_M -= microbe.vf_phag_excreted 
                    self.number_vir_interactions += 1
                else:
                    self.phag_prob_M = self.phag_prob_M_min
                    self.number_vir_interactions += 1

    def die_divide(self):
        if self.alive:
        # UPDATE K_NET
            self.set_netgrowth([x[self.pos[0]][self.pos[1]] for x in self.model.ab_conc_grid_alldrugs])           
        #  decide whether DIVIDE OR DIE based on if K_NET is positive or negative
            if self.K_NET < 0:
                if self.model.random.random() < abs(self.model.stepsize*self.K_NET):
                    self.die()
            elif self.K_NET > 0:
                if self.model.random.random() < (self.model.stepsize*self.K_NET):
                    self.divide()
      
    def die(self):
        self.model.kill_microbes_list.append(self) # add to kill list
        self.model.number_microbes -= 1
        self.model.number_microbes_species[self.species_name] -= 1
#        self.model.number_per_species[[i for i in range(self.model.number_species) if self.species_name==self.model.species_types[i]][0]] -= 1        
        self.alive = False
     
    def divide(self):
        # reuse the set agent characteristics back to imitate newborn agent with or without inheriting some of the charactersistics, 1st daughter
        self.reset_microbe() # reset and move to random patch
        # create one more new agent of the same species as copy of 1st daughter
        self.model.create_copy(self)
        
    def reset_microbe(self):
        # Parameters specific to each agent
        # self.microbe_id does not change
        # reset interaction effects
        self.MIC = copy.deepcopy(self.MIC_original)# deepcopy list! # current MIC of this agent, initialized to be species specific MIC
        self.k_repmax = self.k_repmax_original

    def move_rand(self):
        coords = (self.model.random.randrange(0, self.model.gridrange), self.model.random.randrange(0, self.model.gridrange))
        self.model.grid.move_agent(self, coords)

    # move self to position pos        
    def move_daughter(self,pos):
#        self.model.number_microbes_grid[self.pos[0]][self.pos[1]] -=1 # remove counter from old position
        self.model.grid.move_agent(self,(pos[0],pos[1]))
#        self.model.number_microbes_grid[self.pos[0]][self.pos[1]] +=1 # add counter at new position        


##################################################################################################################

        
class MyNeutrophil(Agent): # neutrophils

# AGENT INITIALIZATION   
    
    def __init__(self, model, neutrophil_id, neutrophils):
        super().__init__(neutrophil_id, model)
        
        # neutrophil specific parameters, these parameters do not change during the simulation
        self.neutrophil_id = neutrophil_id # unique agent id
        self.phag_capacity = 48 # total number of bacteria neutrophil can phagocytose # vandebroucke-grauls et al., 1984
        self.phag_capacity_per_min = 6 # number of bacteria neutrophil can phagocytose per min # vandebroucke-grauls et al., 1984
        self.half_life = neutrophils[1]        
        
        # parameters specific to each neutrophil
        self.alive = True 
        self.recruit_time = self.model.time # time of recruitment of the neutrophil
        self.death_time = np.nan
        self.phag_prob_N = neutrophils[2] * self.model.stepsize # first order phagocytose rate constant * delta t
        self.age = self.model.random.random() # neutrophils dies when age >= 1, updated every step
        self.dying_age = self.model.random.random()    

        # counters
        self.number_phagocytosed = 0 # counter for total number of microbes phagocytosed
        self.phagocytosed_per_step = 0 # counter for number of microbes phagocytosed, resets every step
   
# STEP FUNCTION FOR NEUTROPHILS 
        
    def step(self):
        if self.alive == True:
            self.move()
            self.lifeupdate()
  
# ACTION FUNCTIONS FOR NEUTROPHILS        

    def move(self): # move to one of 8 neighbouring patches or stay on the same patch, all with equal chances
        possible_steps = self.model.grid.get_neighborhood(
            self.pos,
            moore=True,
            include_center=True)
        new_position = self.random.choice(possible_steps)
        self.model.grid.move_agent(self, new_position)  
        self.stepid_last_move = self.model.step_id
        
    def lifeupdate(self):
        # update age
        if self.model.time != self.recruit_time:    
            self.age += self.model.stepsize 

        # phagocytose?
        self.phagocytosed_per_step = 0 # reset counter for number microbes phagocytosed per step  
        
        microbe_cellmates = self.model.give_microbe_cellmates([self.pos]) # gives a list of all microbes in current patch
        self.random.shuffle(microbe_cellmates) # shuffles the list of microbes so there is no preference for microbes with lowest id
        
        for microbe in microbe_cellmates: 
            if microbe.alive and self.can_phagocytose(microbe):
                self.phagocytose(microbe)                        
               
        # die?
        if np.e**-(self.age*(np.log(2)/self.half_life)) < self.dying_age or self.number_phagocytosed >= self.phag_capacity: 
            self.die()

    def can_phagocytose(self, microbe): # returns True if microbe and neutrophil meet all conditions for phagocytosis
        phag_capacities = self.phagocytosed_per_step < (self.phag_capacity_per_min*self.model.stepsize) and self.number_phagocytosed < self.phag_capacity # if either of the phagocytose capacities is reached this returns False
        phag_prob_M =  microbe.phag_prob_M > self.model.random.random() # if the probability of phagocytosis is high, higher chance this returns True
        phag_prob_N = self.phag_prob_N > self.model.random.random() # if the probability of phagocytosis is high, higher chance this returns True
        return phag_capacities and phag_prob_M and phag_prob_N
                    
    def phagocytose(self, microbe): # phagocytosis of microbe by neutrophil
        microbe.die()
        self.number_phagocytosed += 1
        self.phagocytosed_per_step += 1
        #print("- neutrophil {} killed microbe {}".format(self.neutrophil_id, microbe.microbe_id))
        
    def die(self): # death of neutrophil
        self.model.kill_neutrophils_list.append(self) # add neutrophil to model kill list 
        self.model.number_neutrophils -= 1
        self.death_time = self.model.time
        self.alive = False

        
##################################################################################################################
        

class MyModel(Model):

# MODEL INITIALIZATION

    def __init__(self, microbe_n, microbe_capacity, all_species, all_species_virulence, # microbes
                 neutrophil_n, neutrophil_capacity, neutrophils, # neutrophils
                 ab_conc_time, drugs, runtime, gridrange, # antibiotics
                 move_type="stepwise", n_steps = 1, capacity_type = 'global', replication_radius = 0, data_collection = 10, data_agents = False, activation = "random", stepsize = 1):
        super().__init__()
        self.running = True # necessary for batch run, is set to False when all bacteria have died, This stops the simulation run

        # simulation parameters
        self.runtime = runtime
        self.stepsize = stepsize # in mins
        self.time = 0
        self.step_id = 0
        self.data_collection = data_collection # interval for datacollection in mins
        self.data_agents = data_agents # input parameter that decides whether data from individual agents is saved or not
        self.starting_number_microbes = microbe_n
        
        if activation == "random":
            self.schedule = RandomActivation(self)
        if activation == "staged": # USE RANDOM SETTING ONLY
            self.stage_list = ["move_agent", "interact", "divide_die", "update_interaction_effects"]
            self.schedule = StagedActivation(self, self.stage_list, shuffle = True, shuffle_between_stages = False)

        # model parameters
        self.microbe_capacity = microbe_capacity # maximal number of microbe agents
        self.capacity_type = capacity_type # USE capacity_type='global'
        self.replication_radius = replication_radius # irrelevant if capacity_type='global'        

        # spatial grid and antibiotic concentration
        self.gridrange = gridrange
        self.grid = MultiGrid(self.gridrange, self.gridrange, torus=True)
        self.move_type = move_type # input parameter that decides if agent movement is stepwise (one patch at a time) or random
        self.n_steps = n_steps # number of move steps each agent does per timestep, default=1    
        
        #create grid for initial distribution of antibiotic concentration
        self.ab_conc_time = ab_conc_time # applied concentration over time as predefined, not on grid yet
        self.set_ab_conc()
        
        # define drug characteristics
        self.G_min=drugs[0] # G_min: minimal net growth rate reached, can be negative for kill
        self.k = drugs[1] # Hill factor
        self.drugtype = drugs[2] # cidal or static
        
        # unique microbe names and some running counters
        self.number_microbes = 0 # counts alive agents        
        self.microbe_id = 0.1 # to assign unique name, .1 stands for microbe agent, .2 stand for neutrophil agent
        self.number_species = len(all_species) # number of species        
        self.species_names = [species[0] for species in all_species]
        self.number_microbes_species = dict.fromkeys([species[0] for species in all_species], 0) # dict of species names to count number of alive agents for every species
#        self.number_per_species = [0 for x in range(self.number_species)] # counts agents per species         
        self.species_types = [x[0] for x in all_species] # list of all the species types present
        self.nm = copy.deepcopy(self.number_microbes_species)          
#        self.number_microbes_grid = [[0 for j in range(self.gridrange)] for i in range(self.gridrange)] # counts agents per patch

        # neutrophil parameters and counters
        self.neutrophil_capacity = neutrophil_capacity # maximum number of neutrophil agents
        self.neutrophils = neutrophils # max_recruit_rate, halflife, phag_prob_N
        self.max_recruit_rate = neutrophils[0] # max recruitment rate of neutrophils, per min 
        self.number_neutrophils = 0 # counts alive neutrophils        
        self.number_recruited = 0 # counts how many neutrophils are recruited
        
        # model level virulence
        self.vf_recruit_species = dict(zip((self.species_names),(1-all_species_virulence[i][2] for i in range(len(all_species_virulence))))) # dict of species names and recruitment virulence effects (vf_recruit)
        
        # initialize
        self.kill_microbes_list = list() # initialize microbe kill list    
        self.kill_neutrophils_list = list() # initialize neutrophil kill list
        
    # CREATE MICROBES
        
        # make a random starting ditribution of species with approximately equal numbers each, ratio is always 1:1 for all species
        x = self.starting_number_microbes/self.number_species
        ns = [self.random.gauss(x,x*0.1) for i in range(self.number_species)] # draw random numbers between 0 and 1, uniform disctribution
        ns = [ i/sum(ns) for i in ns ] # divide each element by sum
        ns = [int(i*self.starting_number_microbes) for i in ns] # multiply each element by desired sum, this leads to approximately n_init bacteria at the start
        ns = [x+1 if x==0 else x for x in ns]
            
        # create agents
        for i in range(self.number_species):
            self.create_microbes(ns[i], all_species[i], all_species_virulence[i])#, [x[i] for x in all_interactions], [x[i] for x in all_interactions_env]) # create chosen number of each type of species

    # CREATE NEUTROPHILS

        self.neutrophil_id = 0.2 # to assign unique name, .1 stands for microbe agent, .2 stand for neutrophil agent
        self.starting_number_neutrophils = neutrophil_n
        self.create_neutrophils(self.starting_number_neutrophils)

    # DATA COLLECTION
        
        model_reporters = {
            "timestep": "step_id", 
            "Time": "time", 
            "Number_microbes": "number_microbes",
            "Number_species": "nm",
            "Number_neutrophils": "number_neutrophils",
            "Number_recruited": "number_recruited",
            "AB_fix":"ab_conc"
            }
        agent_reporters = {
             # Microbe reporters
             "Species": lambda a: a.species_name if hasattr(a, "species_name") else None, 
             "Growthrate": lambda a: a.replicationrate if hasattr(a, "replicationrate") else None,            
             "Killrate": lambda a: a.deathrate if hasattr(a, "deathrate") else None,
             "NetGrowthrate": lambda a: a.K_NET if hasattr(a, "K_NET") else None,
             "Phag_prob_M": lambda a: a.phag_prob_M if hasattr(a, "phag_prob_M") else None,
             "Number_vir_interactions": lambda a: a.number_vir_interactions if hasattr(a, "number_vir_interactions") else None,
             
             # Neutrophil reporters
             "Recruit_time": lambda a: a.recruit_time if hasattr(a, "recruit_time") else None,
             "Age": lambda a: a.age if hasattr(a, "age") else None,
             "Death_time": lambda a: a.death_time if hasattr(a, "death_time") else None,
             "Number_phagocytosed": lambda a: a.number_phagocytosed if hasattr(a, "number_phagocytosed") else None,
             "Phag_prob_N": lambda a: a.phag_prob_N if hasattr(a, "phag_prob_N") else None
             }    

        # set up data collection during runtime 
        if self.data_agents: #only collect agent reporters if data_agents is True
            self.datacollector = DataCollector(model_reporters = model_reporters, agent_reporters = agent_reporters)
        else: 
            self.datacollector = DataCollector(model_reporters = model_reporters)
            
        self.collect_data() # collect initial data before simulation starts, for time=0    
 
# FUNCTIONS FOR MODEL INITIALIZATION
        
    def set_ab_conc(self):
        self.ab_conc = [x[self.step_id] for x in self.ab_conc_time]
        self.ab_conc_grid_alldrugs=[[[self.ab_conc[drug_number] for j in range(self.gridrange)] for i in range(self.gridrange)] for drug_number in range(len(self.ab_conc))]

# MODEL FUNCTIONS
    
    def step(self):

        # UPDATE time and conc for next step
        self.step_id += 1
        self.time += self.stepsize
        self.set_ab_conc() # set antibiotic concentration to next concentration and reset ab_change to zero    
        
        # 1) STEP AGENTS
        self.schedule.step() # activate each agents in random order and then run the agent's step function
        
        # 2) EXECUTE KILLING
        self.kill_microbes() # execute to kill all the agents that have been placed on the kill list
        # this has to be done on the model level and cannot be done on the agent level
        
        # 3) STOP?
        # stop simulation run when all bacteria have died
        if self.number_microbes == 0:
            self.running = False
            
        # 4) RECRUIT NEUTROPHILS 
        # neutrophils are recruited after every model STEP AGENTS
        # recruit neutrophils based on number_microbes and number_neutrophils
        if self.running == True: 
            self.recruit_neutrophils()          

        # 4) COLLECT DATA 
        # save variables for each time step for data collection
        # collect data approximately every x minutes as defined in data_collection
        # collect also for the last time step, both when stopped because no cells are left or when max_steps is reached
        if self.time%self.data_collection < self.stepsize or self.running==False or round(self.runtime/self.stepsize)==self.step_id: 
            self.collect_data()
            self.kill_neutrophils() # only killed after data collection points because dead neutrophil agent still carry data of interest
            #print(self.time, self.number_microbes, self.number_neutrophils, self.number_recruited, self.number_death_neutrophils)
        
        # 6) RESET NEUTROPHIL COUNTERS
        self.number_recruited = 0 # reset
  
   
    def collect_data(self):
        self.nm=copy.deepcopy(self.number_microbes_species)
        self.datacollector.collect(self)
        
# FUNCTIONS FOR MICROBES

    def create_microbes(self, n, species, virulence_factors):#, interactions, interactions_env): # use onlu during runtime through create_copy to set agent characteristics correctly
        for i in range(n):
            a = MyMicrobe(self, self.microbe_id + 1, species, virulence_factors)#, interactions, interactions_env) # create new agent with unique name
            self.schedule.add(a)
            coords = (self.random.randrange(0, self.gridrange), self.random.randrange(0, self.gridrange))
            self.grid.place_agent(a, coords)
            self.microbe_id += 1
            self.number_microbes += 1
            self.number_microbes_species[species[0]] += 1
#            Cathi way            self.number_per_species[[i for i in range(self.number_species) if species[0]==self.species_types[i]][0]] += 1
        return a    

    def create_copy(self, a): # creates a copy of 1st daughter,  no need to deal with inheritance, only redraw distributions basically
        a_copy = self.create_microbes(1, a.species, a.virulence_factors)#, a.interactions, a.interactions_env) # position is random

        # Species specific parameters, these parameters do not change during the simulation
        a_copy.species = a.species # variable containing most important species characteristics
        a_copy.species_name = a.species_name # name of species
        a_copy.MIC_original = copy.deepcopy(a.MIC_original) # original MIC for this species
        a_copy.k_repmax_original = a.k_repmax_original # maximal growth rate for this species
        
        # Parameters specific to each agent
        a_copy.MIC = copy.deepcopy(a.MIC) # current MIC of this agent, initialized to be species specific MIC
        a_copy.k_repmax = a_copy.k_repmax_original
        a_copy.move_daughter(a.pos)        
    
    def kill_microbes(self): 
        for x in self.kill_microbes_list:
            self.grid.remove_agent(x)
            self.schedule.remove(x)
        self.kill_microbes_list = list() # reset list for next step    
        
# FUNCTIONS FOR NEUTROPHILS

    def create_neutrophils(self, n): 
        if n > 0:
            for i in range(n):
                a = MyNeutrophil(self, self.neutrophil_id, self.neutrophils) # create new agent with unique id
                self.schedule.add(a)
                coords = (self.random.randrange(0, self.gridrange), self.random.randrange(0, self.gridrange))
                self.grid.place_agent(a, coords)
                self.neutrophil_id += 1
                self.number_neutrophils += 1
            return a

    def recruit_neutrophils(self):
        # recruitment v2.2 - recruitment based on only paper of Ternent et al., 2015
        self.neutrophil_restriction = (1 - self.number_neutrophils/self.neutrophil_capacity) # number 0-1, more neutrophils -> less recruitment 
        
        self.virulence_recruit = 0
        for species_name in self.species_names:
            self.virulence_recruit += self.vf_recruit_species[species_name] * (self.number_microbes_species[species_name]/self.number_microbes)

        n = self.max_recruit_rate * self.number_microbes * self.neutrophil_restriction * self.stepsize * self.virulence_recruit
        
        recruit_n = int(n + self.random.random())
        if self.number_neutrophils + recruit_n > self.neutrophil_capacity:
            recruit_n = self.neutrophil_capacity - self.number_neutrophils
        self.create_neutrophils(recruit_n)
        self.number_recruited += recruit_n
    
    def kill_neutrophils(self): 
        for x in self.kill_neutrophils_list:
            self.grid.remove_agent(x)
            self.schedule.remove(x)
        self.kill_neutrophils_list = list() # reset list for next step   
        
# GENERAL AGENT FUNCTIONS

    def give_microbe_cellmates(self, pos):
        microbe_cellmates = []
        cellmates = self.grid.get_cell_list_contents(pos) # get list of present agents in current pos
    
        if len(cellmates) != 0:
            for cellmate in cellmates:
                if isinstance(cellmate, MyMicrobe): # check if cellmates are microbes
                    microbe_cellmates.append(cellmate)
        return microbe_cellmates

    def give_neutrophil_cellmates(self, pos):
        neutrophil_cellmates = []
        cellmates = self.grid.get_cell_list_contents(pos) # get list of present agents in current pos
    
        if len(cellmates) != 0:
            for cellmate in cellmates:
                if isinstance(cellmate, MyNeutrophil): # check if cellmates are neutrophils
                    neutrophil_cellmates.append(cellmate)
        return neutrophil_cellmates  