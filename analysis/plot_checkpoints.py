import pickle
import os 
import sys
import numpy as np 
import matplotlib.pyplot as plt
sys.path.append(os.path.dirname(os.path.abspath('')))
# import importlib 
# import afpo
# importlib.reload(afpo)
from afpo import AFPO


seed = 0
gen = 75


checkpoint_file = f'../checkpoints/run{seed}_gen{gen}.p'
with open(checkpoint_file, 'rb') as f:
    afpo, np_rng_state = pickle.load(f)
    plt.plot(np.max(afpo.fitness_data[:gen,:, 0], axis=1))
# print(afpo.return_best().voxel_id_1, afpo.return_best().voxel_id_2, afpo.return_best().voxel_id_3, afpo.return_best().voxel_id_4)
plt.xlabel('Generation')
plt.ylabel('Shear Modulus')
plt.grid()
plt.savefig("fitness.png")