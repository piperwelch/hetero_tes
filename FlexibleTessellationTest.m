%% Flexible Tessellation
%{
NOTE: This is all the code which has changed to add more control to tessellation 
patterns. This code and Jamming.m was the only code edited.   
%}
%%
het_dir = 'Heterogenous/';

% read files of single configurations: 
dir_str = strcat(het_dir, 'Data/theta0.7854');
pos_unique = readmatrix([dir_str,'/pos_N6.txt']);
Gn = readmatrix([dir_str,'/Gn_N6.txt']);
N = 6;

% initialize variables & arrays
Kp = 1;
Kw_ratio = 0.1;
Kpw_ratio = 1;
Kw_ratio_tess = Kw_ratio * 0.5; %this is currently half to agree with single voxel shear modulus
%%
tic
cd([het_dir 'Jamming']);
Nvoxel_row = 4;
Nvoxel = Nvoxel_row * Nvoxel_row;
xy_voxels = zeros(size(pos_unique, 1), 2*N+8);
D_every = pos_unique(:, 1);

%rotations single voxels to all have bottom edge as horizontal (mostly formality)
for j = 1:size(pos_unique, 1)
        xy_voxels(j, :) = RotateConfig_Voxel(N, pos_unique(j, 2:2*N+9), 1, 4,  [1; 2; 3; 4]);
end


% create conid array - elements are index of config you want in that spot in array
tessname = '11_config1';
conid = ones(Nvoxel, 1);
conid(1) = 1; 
conid(2) = 2; conid(3) = 5; conid(4) = 4; conid(5) = 5; conid(6) = 6; 
conid(7) = 7; conid(8) = 8; conid(9) = 9; conid(10) = 10; conid(11) = 11; 
conid(12) = 12; conid(13) = 13; conid(14) = 14; 
% tessellate system -> input conid, jamming uses it to place single voxels into array
test = 0; L0 = 1; 
%disp(Kw_ratio)
[Wlist, linklist, linklist_inner, ext_list] = ConstructTruss(Nvoxel_row);
[xy_p_all, xy_c_all, L0_voxel, D_all] = Jamming(conid, xy_voxels, N, D_every, Gn, L0, Nvoxel_row, Kp, Kw_ratio_tess, Kpw_ratio, test);

% plot configuration to check: 
PlotConfig(xy_p_all, xy_c_all, N, Nvoxel_row, D_all, Wlist, linklist, ext_list, 1)
toc
%%
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

Ntheta = 30;
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

%plot(theta_all / pi, Gn(1, :), 'x-') % single voxel G indexed by unique config
%plot(theta_all / pi, Gn(2, :), 'x-')
disp(Gn)
xlabel("\theta / \pi")
ylabel("G_{tess}")

legend("tess")

%% extract G_amp, G_dc, theta_0 for tessellated

tess_info = FitSines(G_step', theta_all', 1); % note: inputs need to be N_trajectories x N_angles
tess_Gamp = tess_info(4); % G_amp; as in G = G_amp * sin(4*th + th_0) + G_dc
tess_Gdc = tess_info(5); % G_dc
tess_th0 = tess_info(6); % theta_0 (in radians)
disp(tess_Gdc)