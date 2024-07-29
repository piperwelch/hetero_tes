function [sigma_shear, gamma_all, G_all, xy_p_t, xy_c_t] = ShearModulus(xy_p_all, xy_c_all, N, Nvoxel_row, D_all, L0, Kp, Kpw_ratio, Kw_ratio, dgamma, theta)
%% construct wall vertex list
Nvoxel = Nvoxel_row * Nvoxel_row;
Nvert = (Nvoxel_row + 1)^2;
Np = Nvoxel * N;
[Wlist, linklist, linklist_inner, ext_list] = ConstructTruss(Nvoxel_row);
%% rotate configuration
[xy_p_all, xy_c_all] = RotateConfig(xy_p_all, xy_c_all, Np, Nvert, theta);
xy_c_ori = reshape(xy_c_all, Nvert, 2)';
%% FIRE parameters
%Kp = 1;
Kpw = Kpw_ratio * Kp;
Kw = Kw_ratio * Kp;
% Kpw = Kpw_ratio * Kw;
% Kp = Kp_ratio * Kw;

dt_fire = 0.01 / sqrt(max([1, 4 * Kpw, Kw]));
Fthresh = 1E-14 * max([1, 4 * Kpw, Kw]);
Nt_fire = 1E7;
%%
Nstep = 20;
xy_p_t = zeros(2 * Np, Nstep + 1);
xy_c_t = zeros(2 * Nvert, Nstep + 1);
xy_p_t(:, 1) = xy_p_all;
xy_c_t(:, 1) = xy_c_all;

gamma_all = (0:Nstep)' * dgamma;
sigma_shear = zeros(Nstep + 1, 4);
sigma_shear_pp = zeros(Nstep + 1, 4);
sigma_shear_pw = zeros(Nstep + 1, 4);
sigma_shear_ww = zeros(Nstep + 1, 4);

[sigma, sigma_pp, sigma_pw, sigma_ww, ~] = StressTensor([xy_p_all; xy_c_all], N, Nvoxel, Nvert, D_all, L0, Wlist, linklist, Kp, Kpw, Kw);
sigma_shear(1, :) = reshape(sigma, 1, 4);
sigma_shear_pp(1, :) = reshape(sigma_pp, 1, 4);
sigma_shear_pw(1, :) = reshape(sigma_pw, 1, 4);
sigma_shear_ww(1, :) = reshape(sigma_ww, 1, 4);

for it = 1:Nstep
    %fprintf('it: %d\n', it);
    gamma_it = it * dgamma / 2;
    J = [1, gamma_it;
         gamma_it, 1];

    xy_c_strain = J * xy_c_ori;
    xy_c_all = reshape(xy_c_strain', 2 * Nvert, 1);
    [xy_p_all, xy_c_all, ~] = FIRE_Disk(xy_p_all, xy_c_all, N, Nvoxel, Nvert, D_all, L0,...
                                        Wlist, linklist, ext_list, Kp, Kpw, Kw, Fthresh, dt_fire, Nt_fire);
    xy_p_t(:, it + 1) = xy_p_all;
    xy_c_t(:, it + 1) = xy_c_all;

    [sigma, sigma_pp, sigma_pw, sigma_ww, ~] = StressTensor([xy_p_all; xy_c_all], N, Nvoxel, Nvert, D_all, L0, Wlist, linklist, Kp, Kpw, Kw);
    sigma_shear(it + 1, :) = reshape(sigma, 1, 4);
    sigma_shear_pp(it + 1, :) = reshape(sigma_pp, 1, 4);
    sigma_shear_pw(it + 1, :) = reshape(sigma_pw, 1, 4);
    sigma_shear_ww(it + 1, :) = reshape(sigma_ww, 1, 4);
end

G_fit = polyfit(gamma_all, sigma_shear(:, 2), 1);
G_fit_pp = polyfit(gamma_all, sigma_shear_pp(:, 2), 1);
G_fit_pw = polyfit(gamma_all, sigma_shear_pw(:, 2), 1);
G_fit_ww = polyfit(gamma_all, sigma_shear_ww(:, 2), 1);
G_all = -[G_fit(1); G_fit_pp(1); G_fit_pw(1); G_fit_ww(1)];