

class Material: 
    def __init__(self, random, id):
        # self.random = random
        self.id = id
        self.age = 0 
        self.make_genome(random)


    def make_genome(self, random):
        num_voxels = 15
        
        self.genome = random.randint(1, num_voxels, 16)


    def mutate(self, random):
        num_configs = 15
        mutation_index = random.randint(0, 15) #choose voxel 1-16 to mutate 
        new_voxel = random.randint(1, num_configs) #choose new voxel from 1-14
        
        new_genome = self.genome.copy()
        new_genome[mutation_index] = new_voxel
        
        self.genome = new_genome.copy()


    def print(self, verbose):
        print(self.fitness, self.genome)
            