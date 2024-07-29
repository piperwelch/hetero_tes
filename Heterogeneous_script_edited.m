%%
current_dir = pwd;
single_voxel_dir = [current_dir '/Heterogenous/Single/'];
het_dir = [current_dir '/Heterogenous/'];
packing_dir = [current_dir '/Heterogenous/UniquePacking/'];
data_dir = [current_dir '/Heterogenous/Data/'];

%% generate many single voxel configs
cd([single_voxel_dir 'Jamming']);

L0 = 1;
N = 4;
N_configs = 100;
pos_all = zeros(N_configs, 1 + 2 * N + 8);
xy_voxels_old = zeros(2, 2 * N + 8);
D_every = zeros(N_configs, 4);

%interaction energy definitions - based on Zhang et. al. paper numbers 
Kp = 1;
Kw_ratio = 0.1;
Kpw_ratio = 1;


%generate N_configs random voxels with fixed boundary conditions (grows particles to jam)    
tic
for i = 1:N_configs
    fprintf("on config %.0f \n", i)
    N = 4; pid = i; ratio = 1.4; Dtype = 1; test = 0; theta = pi/3;
    [xy_voxel, D_all, D] = Jamming(N, pid, ratio, Dtype, Kp, Kw_ratio, Kpw_ratio, test, theta);
    pos_all(i, 1) = D; % diameter
    D_every(i, :) = D_all;
    pos_all(i, 2:2*N+9) = xy_voxel;
end
toc
% takes roughly ~5-10 sec per config

%% 

% save('100_60deg_configs.mat', 'pos_all', 'D_every', 'Kp', 'Kw_ratio', 'Kpw_ratio')
load([data_dir '100_90deg_configs.mat']) % load in previously saved data

%% sample plot call

cd([single_voxel_dir 'Jamming']);
p_ind = 8;
PlotConfig(pos_all(p_ind, 2:2*N+9), N, D_every(p_ind, :), 1)

%% unique configurations - picks out unique configurations from random voxels (doesn't currently sort out rattlers)

cd(packing_dir);
[pos_unique_act, pos_unique, config_idx, appear] = ObtainDistinctPacking(N, pos_all);

% write unique config info to text file
writematrix(pos_unique_act, [data_dir 'unique_pos_90deg.txt'])
writematrix(config_idx, [data_dir 'idx_90deg.txt'])

%% plots all unique configurations found in prev section (new figs)

cd([single_voxel_dir 'Jamming']);
Nh = floor(N / 2);
Gn = [ratio * ones(Nh, 1); ones(Nh, 1)];

% for i = 1:max(config_idx)
%     %fprintf('i = %d\n', i)
%     xy_voxel = pos_unique(i, 2:2*N + 9)';
%     D_all = pos_unique(i, 1) * Gn;
%     PlotConfig(xy_voxel, N, D_all, 1)
%     title(['i = ' num2str(i)]) 
% end

%% single voxels FLW shear modulus
cd([single_voxel_dir 'ShearModulus_Out']);

% ~1sec to run 50 thetas for 2 configs, single vox
step = 50;  % number of angles
dgamma = 5e-9;

N_unique = size(pos_unique_act, 1);
theta_all = linspace(0, pi / 2, step)';
full_G = zeros(N_unique, step); % unique G vs thetas;
G_step = zeros(step, 1);

color = ['o'; '*'; 'x'; 'v'; 's'; 'p'; 'o'];

figure
tic
for j = [1:N_unique]
    for i = 1:step
%         [sigma_all_shear, sigma_p_shear, sigma_w_shear, Fext_all, gamma_all, G_all, xy_t] = ...
%                 ShearModulus_Strain(xy_voxels(j,:)', L0, N, D_every(j, :), Kp, Kw_ratio, Kpw_ratio, dgamma, theta_all(i));
        try
        unique_idx = unique(config_idx);
        class_inds = find(unique_idx(j)==config_idx);
        [sigma_all_shear, sigma_p_shear, sigma_w_shear, Fext_all, gamma_all, G_all, xy_t] = ...
                ShearModulus_Strain(pos_unique_act(j,2:end)', L0, N, D_every(class_inds(1), :), Kp, Kw_ratio, Kpw_ratio, dgamma, theta_all(i));        
        catch
            continue
        end
        G_step(i) = G_all(1);
    end
    
    full_G(j, :) = G_step;
    
%     title('Single FLW Voxels v. Grid Tessellation, grid \epsilon_{pp} = 0.5')
%     plot(theta_all / pi, G_step, 'Marker',color(j), 'DisplayName',['configuration' num2str(j)])
%     xlabel('\theta/\pi')
%     ylabel('G')
%     hold on 
end
toc

%% visualize sanity check
figure(1), clf, hold on

for jj = [1:N_unique]
% for jj = [1:5]
    plot(theta_all, full_G(jj, :), 'o-')
end

%% write shear moduli (corresponding to unique configs) to text file

writematrix(full_G, [data_dir 'unique_G_90deg.txt'])

%% tessellation jamming - uses flexible boundary conditions (shrinks box to jam) 

tic
cd([het_dir 'Jamming']);
Nvoxel_row = 4;
N_grid = 2; %number of different configs in checkerboard
xy_voxels = zeros(N_grid, 2*N+8);

%initialize 2 voxels in checkerboard  
xy_voxels_old(1, :) = pos_unique_act(3, 2:2*N+9);
xy_voxels_old(2, :) = pos_unique_act(2, 2:2*N+9);

D_every(1, :) = pos_unique_act(3, 1) * Gn;
D_every(2, :) = pos_unique_act(2, 1) * Gn;
 
for i = 1:N_grid
    xy_voxels(i, :) = RotateConfig_Voxel(N, xy_voxels_old(i, :), 1, 4,  [1; 2; 3; 4]);
end

test = 0; L0 = 1;
[xy_p_all, xy_c_all, L0_voxel, D_all] = Jamming(N_grid, xy_voxels, N, D_every, L0, Nvoxel_row, Kp, Kw_ratio, Kpw_ratio, test);
[Wlist, linklist, linklist_inner, ext_list] = ConstructTruss(Nvoxel_row);
PlotConfig(xy_p_all, xy_c_all, N, Nvoxel_row, D_all, Wlist, linklist, ext_list, 1)
toc

% really fast, ~1 sec (system is already mostly energy minimized)

%% tessellation compress
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

%% run g vs p || full tess

tic

cd([het_dir 'ShearModulus']);
stepnum = size(P_comp, 1);
G_step = zeros(stepnum, 4);

dgamma = 5e-9;
theta = 0;

%Kpw_ratio = 0.25; Kw_ratio = 0.1;
for i = 1:stepnum
    [sigma_shear, gamma_shear, G_all, xy_p_shear, xy_c_shear] = ...
            ShearModulus(xy_p_all_comp(:, i), xy_c_all_comp(:, i), N, Nvoxel_row, D_all, L0, ...
            Kp, Kpw_ratio, Kw_ratio, dgamma, theta);
    G_step(i, :) = G_all';
end
toc

% plot 
figure
plot(P_comp, G_step(:, 1), 'bo', 'markersize', 8, 'linewidth', 1)
% set(gca,'xscale','log')
xlabel('P')
ylabel('G')

% ~ 15 seconds for 50 pressure values

%% single voxels FXW shear modulus 
% cd([single_voxel_dir 'ShearModulus_FXW']);
% 
% step = 101;
% dgamma = 5e-9;
% 
% theta_all = linspace(0, pi / 2, 101)';
% G_step = zeros(101, 1);
% 
% color = ['o'; '*'];
% 
% figure
% 
% for j = 1:N_grid
%     for i = 1:101
%         [sigma_all_shear, sigma_p_shear, Fext_all, gamma_all, G_all, xy_t] = ...
%                 ShearModulus_Strain(xy_voxels(j,:)', L0, N, D_every(j, :), Kp, Kw_ratio, Kpw_ratio, dgamma, theta_all(i));
%         G_step(i) = G_all(1);
%     end
% 
%     title('G of Single FXW Voxels compared to Grid Tessellation')
%     plot(theta_all / pi, G_step, 'Marker',color(j), 'DisplayName',['configuration' num2str(j)])
%     xlabel('\theta/\pi')
%     ylabel('G')
% 
%     hold on 
% end

%% generate G(theta) for entire tessellation
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

plot(theta_all / pi, G_step, 'gx', 'DisplayName', 'tessellated')
xlabel('\theta/\pi')
ylabel('G')
% hold off
toc

legend()

%%

figure(1), clf, hold on
plot(theta_all / pi, G_step, 'x-')

%% extract G_amp, G_dc, theta_0

tess_info = FitSines(G_step', theta_all', 1); % note: inputs need to be N_trajectories x N_angles
tess_Gamp = tess_info(4); % G_amp; as in G = G_amp * sin(4*th + th_0) + G_dc
tess_Gdc = tess_info(5); % G_dc
tess_th0 = tess_info(6); % theta_0 (in radians)
