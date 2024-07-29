import matlab.engine
import numpy as np 
from afpo import AFPO 

seed = 0 

eng = matlab.engine.start_matlab()
current_dir = eng.pwd()
het_dir = eng.set_file_locations(current_dir, nargout=1)

# Load configs
eng.cd(current_dir)
pos_unique, Kp, Kw_ratio, Kpw_ratio, Kw_ratio_tess, N,  Gn = eng.load_configs(het_dir, nargout=7)

print('engine set')
gens = 100
pop_size = 50
afpo = AFPO(seed, gens, pop_size, current_dir, het_dir, pos_unique, Kp, Kw_ratio, Kpw_ratio, Kw_ratio_tess, N,  Gn)
afpo.run(eng)