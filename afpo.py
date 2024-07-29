'''
Created on 2024-06-20 10:40:59
@author: piperwelch 

Description: Implementation of Age-Fitness Pareto Optimization (AFPO) 
https://dl.acm.org/doi/10.1145/1830483.1830584

NOTE: For AFPO, population size must be significantly large so that not all individuals in the population are on the pareto front. 
'''

import numpy as np
import copy
import operator
import pickle
import os
# import random
from material import Material 
# import sys
import numpy as np
# import matplotlib 
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# from MD_functions import FIRE_VL
# from ShearModulus import ShearModulus
# from PlotFunctions import ConfigPlot_DiffSize
from concurrent.futures import ProcessPoolExecutor


class AFPO:
    def __init__(self, random_seed, gens, pop_size, current_dir, het_dir, pos_unique, Kp, Kw_ratio, Kpw_ratio, Kw_ratio_tess, N,  Gn, checkpoint_every=1):
        self.Kp = 1.0  
        self.Kw_ratio = 0.1
        self.Kpw_ratio = 1.0  
        self.test = 0 
        self.seed = random_seed
        # eng = eng
        # Seed rng 
        np.random.seed(self.seed)
        # random.seed(self.seed)
        self.current_dir = current_dir
        self.gens = gens
        self.pop_size = pop_size
        self.checkpoint_every = checkpoint_every
        self.next_id = 0
        self.fitness_data = np.zeros(shape=(self.gens+1, self.pop_size, 2))
        self.het_dir = het_dir
        self.N = N
        self.Gn = Gn
        self.pos_unique = pos_unique
        self.Kp = Kp
        self.Kw_ratio = Kw_ratio
        self.Kpw_ratio = Kpw_ratio
        self.Kw_ratio_tess = Kw_ratio_tess
        os.makedirs('checkpoints/', exist_ok=True)
        self.create_initial_population()

    
    def create_initial_population(self):
        self.pop = [] # population is a list of materials
        for i in range(self.pop_size): # initialize random materials
            self.pop.append(Material(random=np.random, id=self.next_id))
            self.next_id+=1 


    def run(self, eng, continue_from_checkpoint=False, additional_gens=0):
        if continue_from_checkpoint:
            max_gens = self.curr_gen + additional_gens

            # Expand fitness data matrix to account for additional gens
            new_fitness_data = np.zeros((max_gens + 1, self.pop_size, 2))

            new_fitness_data[0:self.curr_gen,:,:] = self.fitness_data[0:self.curr_gen,:,:]
            self.fitness_data = new_fitness_data

            for i,material in enumerate(self.pop):   
                self.fitness_data[self.curr_gen,i,0] = material.fitness
                self.fitness_data[self.curr_gen,i,1] = material.age
            self.curr_gen += 1

            self.gens = max_gens
            while self.curr_gen < self.gens + 1:

                self.perform_one_generation(eng)
                if self.curr_gen % self.checkpoint_every == 0:
                    self.save_checkpoint()

                print("GEN: {}".format(self.curr_gen))
                if self.curr_gen % 5 == 0:
                    self.print_best(verbose=False)
                
                self.curr_gen += 1
        else:
            self.curr_gen = 0
            self.evaluate_generation_zero(eng)
                        
            while self.curr_gen < self.gens + 1: # Evolutionary loop
                
                self.perform_one_generation(eng)
                if self.curr_gen % self.checkpoint_every == 0:
                    self.save_checkpoint()
                
                print("GEN: {}".format(self.curr_gen), flush=True)
                self.print_best(verbose=False)
                self.curr_gen += 1

        return self.return_best(), self.fitness_data


    def evaluate_generation_zero(self, eng):
        
        # Evaluate individuals in the population
        
        self.evaluate(eng, self.pop)

        for i,material in enumerate(self.pop):            
            # Record fitness statistics    
            self.fitness_data[self.curr_gen,i,0] = material.fitness
            self.fitness_data[self.curr_gen,i,1] = material.age

        print("GEN: {}".format(self.curr_gen))
        self.print_best(verbose=False)

        self.curr_gen += 1


    def perform_one_generation(self, eng):

        self.increase_age()
        children = self.breed()
        children = self.insert_random(children)

        # Evaluate children
        self.evaluate(eng, children)
        for child_material in children:
            # Extend population by adding child material (extends to pop_size*2+1 individuals every generation then gets reduced back to pop_size)
            self.pop.append(child_material) 

        self.survivor_selection()

        # Record statistics 
        for i, material in enumerate(self.pop):
            self.fitness_data[self.curr_gen,i, 0] = material.fitness
            self.fitness_data[self.curr_gen,i, 1] = material.age


    def increase_age(self):
        for material in self.pop:
            material.age += 1


    def breed(self):
        children = []
        for i in range(self.pop_size):

            # # Parent Selection via Tournament Selection (based on fitness only)
            parent = self.tournament_selection()
            
            # # Create offspring via mutation
            child = copy.deepcopy(self.pop[parent])
            child.id = self.next_id
            self.next_id += 1
            child.mutate(np.random)
            children.append(child)

        return children


    def insert_random(self, children):
        children.append(Material(random=np.random, id=self.next_id))
        self.next_id += 1
        return children


    def tournament_selection(self):
        p1 = np.random.randint(len(self.pop))
        p2 = np.random.randint(len(self.pop))
        while p1 == p2:
            p2 = np.random.randint(len(self.pop))

        if self.pop[p1].fitness > self.pop[p2].fitness:
            return p1
        else:
            return p2


    def survivor_selection(self):
        # Remove dominated individuals until the target population size is reached
        while len(self.pop) > self.pop_size:

            # Choose two different individuals from the population
            ind1 = np.random.randint(len(self.pop))
            ind2 = np.random.randint(len(self.pop))
            while ind1 == ind2:
                ind2 = np.random.randint(len(self.pop))

            if self.dominates(ind1, ind2):  # ind1 dominates
                
                # remove ind2 from population and shift following individuals up in list
                for i in range(ind2, len(self.pop)-1):
                    self.pop[i] = self.pop[i+1]
                self.pop.pop() # remove last element from list (because it was shifted up)

            elif self.dominates(ind2, ind1):  # ind2 dominates

                # remove ind1 from population and shift following individuals up in list
                for i in range(ind1, len(self.pop)-1):
                    self.pop[i] = self.pop[i+1]
                self.pop.pop() # remove last element from list (because it was shifted up)

        assert len(self.pop) == self.pop_size


    def dominates(self, ind1, ind2):
        # Returns true if ind1 dominates ind2, otherwise false
        if self.pop[ind1].age == self.pop[ind2].age and self.pop[ind1].fitness == self.pop[ind2].fitness:
            return self.pop[ind1].id > self.pop[ind2].id # if equal, return the newer individual

        elif self.pop[ind1].age <= self.pop[ind2].age and self.pop[ind1].fitness >= self.pop[ind2].fitness:
            return True
        else:
            return False
    

    # Parallelize the loop
    def evaluate(self, eng, children):
        for material in children:
            material.fitness = self.evaluate_material(eng, material)
        # self.submit_batch(eng, children)


    def submit_batch(self, eng, children):
        self.children = children 
        
        with ProcessPoolExecutor() as executor:
            results = executor.map(self.evaluate_material, [eng]*len(self.children), self.children)
        print(results)
        for material_index, fitness in results: #todo fitness is difference in G 
            print(material_index, fitness)
            children[material_index].fitness = fitness


    def evaluate_material(self, eng, material):
        
        eng.cd(self.current_dir)
        pos_p_all, pos_c_all, D_all, Wlist, linklist, linklist_inner, ext_list, L0_voxel= eng.jamming(self.het_dir, material.genome, self.pos_unique, self.N, self.Gn, self.Kp, self.Kw_ratio, self.Kpw_ratio, self.Kw_ratio_tess, nargout=8)
        
        #2 seconds

        eng.cd(self.current_dir)
        xy_p_all_comp, xy_c_all_comp = eng.compress(self.het_dir, pos_p_all, pos_c_all, self.N, 4.0, D_all, 1.0, L0_voxel, self.Kp, self.Kpw_ratio, self.Kw_ratio, 0, nargout=2)

        eng.cd(self.current_dir)
        G_step, tess_Gdc = eng.shear_modulus(self.het_dir, xy_p_all_comp,xy_c_all_comp, self.N, 4.0, D_all, 1.0, self.Kp, self.Kpw_ratio, self.Kw_ratio, nargout=2)

        print(tess_Gdc)
        return tess_Gdc


    # def plot_best(self, material):

    #         n_col = 4
    #         n_row = 3
    #         N = (n_col - 1) * n_row + int(np.floor(n_row / 2.0))  # total number of particles

    #         d0 = 1.
    #         y_top_disp = 0 ## compression amount
    #         Lx = d0 * n_col
    #         Ly = (n_row - 1) * np.sqrt(3) / 2 * d0 + d0
    #         Ly -= y_top_disp * d0
    #         xc = np.array([0.0, Lx, Lx, 0.0], dtype = np.float64) # corner x coordinates
    #         yc = np.array([0.0, 0.0, Ly, Ly], dtype = np.float64) # corner y coordinates
    #         k_list = np.ones(N, dtype = np.float64) * 2.

    #         D = np.ones(N, dtype = np.float64) * d0

    #         for particle_index, particle in enumerate(material.particles): #fill in genome 
    #             D[particle_index] = d0 * (1. + particle.expansion) # inflate this specific particle
    #             k_list[particle_index] = particle.stiffness 

    #         x = np.zeros(N, dtype = np.float64)
    #         y = np.zeros(N, dtype = np.float64)

    #         ind = -1
    #         for i_row in range(n_row):
    #             if i_row % 2 == 1:
    #                 n_col_now = n_col
    #             else:
    #                 n_col_now = n_col - 1
    #             for i_col in range(n_col_now):
    #                 ind += 1
    #                 if i_row % 2 == 1:
    #                     x[ind] = (i_col + 0.5) * d0
    #                 else:
    #                     x[ind] = (i_col + 1.) * d0
    #                 y[ind] = i_row * 0.5 * np.sqrt(3) * d0
    #         y = y + 0.5 * d0

    #         mass = np.ones(N, dtype = np.float64)
            
    #         # print(k_list)
    #         FIRE_VL(N, x, y, D, xc, yc, k_list)
    #         ConfigPlot_DiffSize(N, x, y, D, xc, yc, k_list, 1, 0)
    #         # print('hi')
    #         G, strain, stress = ShearModulus(N, x, y, D, xc, yc, k_list)

    #         plt.title(G)
    #         # plt.savefig(f'best_configs_viz/seed{self.seed}_best_material_{material.id}_gen{self.curr_gen}_stiffness{self.change_stiffness}_size{self.change_size}')


    def save_checkpoint(self):

        filename = 'checkpoints/run{}_gen{}.p'.format(self.seed, self.curr_gen)

        # rng_state = random.getstate()
        np_rng_state = np.random.get_state()

        with open(filename, 'wb') as f:
            pickle.dump([self, np_rng_state], f)

    
    def print_population(self, verbose=False):

        for i in range(len(self.pop)):
            self.pop[i].print(verbose=verbose)
        print()


    def print_best(self, verbose):
        best = self.return_best()
        print("BEST material:")
        best.print(verbose=verbose)
        # self.plot_best(best)


    def return_best(self):
        return sorted(self.pop, key=operator.attrgetter('fitness'), reverse=True)[0]