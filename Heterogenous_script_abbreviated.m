%% set file locations
current_dir = pwd;
single_voxel_dir = [pwd '/Heterogenous/Single/'];
het_dir = [pwd '/Heterogenous/'];
packing_dir = [pwd '/Heterogenous/UniquePacking/'];
data_dir = [pwd '/Heterogenous/Data/'];

%% load in old configs, set some parameters manually

% currently have 3 angles of data: {45deg // 60deg // 90deg}
load([data_dir '100_60deg_configs.mat']); % load full jamming data
pos_unique_act = readmatrix([data_dir 'unique_pos_60deg.txt']); % unique configs at theta = 60deg
config_idx = readmatrix([data_dir 'idx_60deg.txt']); % full -> unique indices
full_G = readmatrix([data_dir 'unique_G_60deg.txt']); % unique shear measurements 
N = (size(pos_unique_act, 2)-9)/2;
Nh = floor(N / 2);
ratio = 1.4;
Gn = [ratio * ones(Nh, 1); ones(Nh, 1)];

%% jam checkerboard pattern tessellation

tic
%cd([het_dir 'Jamming']);
Nvoxel_row = 4;
N_grid = 2; %number of different configs in checkerboard
%45 degree 1, 7
%60 degree 4, 5
xy_voxels = zeros(N_grid, 2*N+8);
voxel_1_id = 4; 
voxel_2_id = 5; 
%initialize 2 voxels in checkerboard  
xy_voxels_old(1, :) = pos_unique_act(voxel_1_id, 2:2*N+9); % the first input controls which configs are put into tess.
xy_voxels_old(2, :) = pos_unique_act(voxel_2_id, 2:2*N+9);

D_every(1, :) = pos_unique_act(voxel_1_id, 1) * Gn;
D_every(2, :) = pos_unique_act(voxel_2_id, 1) * Gn;
 
for i = 1:N_grid
    xy_voxels(i, :) = RotateConfig_Voxel(N, xy_voxels_old(i, :), 1, 4,  [1; 2; 3; 4]);
end
disp(xy_voxels)
test = 0; L0 = 1;
[xy_p_all, xy_c_all, L0_voxel, D_all] = Jamming(N_grid, xy_voxels, N, D_every, L0, Nvoxel_row, Kp, Kw_ratio, Kpw_ratio, test);
[Wlist, linklist, linklist_inner, ext_list] = ConstructTruss(Nvoxel_row);
PlotConfig(xy_p_all, xy_c_all, N, Nvoxel_row, D_all, Wlist, linklist, ext_list, 1)
toc

% really fast, ~1 sec (system is already mostly energy minimized)

%% generate tessellation coordinates for different pressures

cd([het_dir 'Compress']);
n_pressure = 40; % number of pressures between 10^-7 and 10^-2 to sample
tic
[xy_p_all_comp, xy_c_all_comp, L0_voxel_comp, P_comp] = Compress(xy_p_all, xy_c_all, N, Nvoxel_row, D_all, L0, L0_voxel, Kp, Kpw_ratio, Kw_ratio, test);
P_t_target = logspace(-7, -2, n_pressure)'; % target pressure values
idx_target = zeros(n_pressure, 1);  % best pressure index
for i = 1:n_pressure
    [~, idx] = min(abs(P_t_target(i) - P_comp));
    idx_target(i) = idx;
end
idx_target = unique(idx_target);
P_comp = P_comp(idx_target);
L0_voxel_comp = L0_voxel_comp(idx_target);
xy_p_all_comp = xy_p_all_comp(:, idx_target');
xy_c_all_comp = xy_c_all_comp(:, idx_target');
toc

% ~ 5 seconds for 50 pressure values

%% generate G(theta) for entire tessellation (at p = 1e-7, can change based on "step" variable)
cd([het_dir 'ShearModulus']);

% ~20 sec to run 50 values of theta
tic
step = 1; % pressure index to measure at
dgamma = 5e-9;

Ntheta = 50;
theta_all = linspace(0, pi / 2, Ntheta)';
G_step = zeros(Ntheta, 1);


for i = 1:Ntheta
    [sigma_shear, gamma_shear, G_all, xy_p_shear, xy_c_shear] = ...
            ShearModulus(xy_p_all_comp(:, step), xy_c_all_comp(:, step), N, Nvoxel_row, D_all, L0, Kp, Kpw_ratio, Kw_ratio, dgamma, theta_all(i));
    G_step(i) = G_all(1);
end

toc


%% visualize G(theta) for tess and single

figure(1), clf, hold on
plot(theta_all / pi, G_step, 'o-')

% plot(theta_all / pi, full_G(1, :), 'x-') % single voxel G indexed by unique config
plot(theta_all / pi, full_G(2, :), 'x-')

xlabel("\theta / \pi")
ylabel("G_{tess}")

legend("tess", "single 1", "single 2")

%% extract G_amp, G_dc, theta_0 for tessellated

tess_info = FitSines(G_step', theta_all', 1); % note: inputs need to be N_trajectories x N_angles
tess_Gamp = tess_info(4); % G_amp; as in G = G_amp * sin(4*th + th_0) + G_dc
tess_Gdc = tess_info(5); % G_dc
tess_th0 = tess_info(6); % theta_0 (in radians)
disp(tess_Gdc)