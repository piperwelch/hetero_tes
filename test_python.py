import matlab.engine
import numpy as np 

np.random.seed(1)
eng = matlab.engine.start_matlab()


current_dir = eng.pwd()
het_dir = eng.set_file_locations(current_dir, nargout=1)


# Load configs
eng.cd(current_dir)
pos_unique, Kp, Kw_ratio, Kpw_ratio, Kw_ratio_tess, N,  Gn = eng.load_configs(het_dir, nargout=7)


eng.cd(current_dir)

# Jam checkerboard pattern tessellation
conid = np.random.randint(1, 15, 16)
print(conid)

pos_p_all, pos_c_all, D_all, Wlist, linklist, linklist_inner, ext_list, L0_voxel= eng.jamming(het_dir, conid, pos_unique, N, Gn, Kp, Kw_ratio, Kpw_ratio, Kw_ratio_tess, nargout=8)

#2 seconds

eng.cd(current_dir)
xy_p_all_comp, xy_c_all_comp = eng.compress(het_dir, pos_p_all, pos_c_all, N, 4.0, D_all, 1.0, L0_voxel, Kp, Kpw_ratio, Kw_ratio, 0, nargout=2)

eng.cd(current_dir)
G_step, tess_Gdc = eng.shear_modulus(het_dir, xy_p_all_comp,xy_c_all_comp, N, 4.0, D_all, 1.0, Kp, Kpw_ratio, Kw_ratio, nargout=2)

print(tess_Gdc)