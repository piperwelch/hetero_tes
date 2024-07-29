function [sigma_all_shear, sigma_p_shear, Fext_all, gamma_all, G_all, xy_t] = ...
                        ShearModulus_Strain(xy_all, L0, N, D_all, Kp, Kpw_ratio, dgamma, theta)
%%
Np_tot = 2 * N;

% xy_all_initial = RotateConfig_Initial(N, xy_all, 1, 2);
xy_all = RotateConfig(N, xy_all, 1, theta);

xy_p = xy_all(1:Np_tot);
xc = xy_all(Np_tot+1:Np_tot+4);
yc = xy_all(Np_tot+5:Np_tot+8);
xy_c = [xc'; yc'];
xy_new = xy_p;

[CL_pp, CL_pw] = ContactList(xy_all, N, D_all);
%%
%Kw = Kp_ratio * Kp;
Kpw = Kpw_ratio * Kp;

dt_fire = 0.01 / sqrt(max([1, 4 * Kpw, 1/D_all(1)^2]));
Fthresh = 1E-14 * max([1, 4 * Kpw, 1/D_all(1)^2]);
Nt_fire = 1E7;

Nstep = 20;
xy_t = zeros(size(xy_all, 1), Nstep + 1);
xy_t(:, 1) = [xy_new; xc; yc];

gamma_all = zeros(Nstep + 1, 1);
sigma_all_shear = zeros(Nstep + 1, 4);
sigma_p_shear = zeros(Nstep + 1, 4);
%sigma_w_shear = zeros(Nstep + 1, 4);
sigma_shear = zeros(Nstep + 1, 1);
Fext_all = zeros(8, Nstep + 1);

[sigma_all, sigma_p] = StressTensorAll_CL(xy_t(:, 1), L0, N, D_all, Kpw, CL_pw);
sigma_all_shear(1, :) = reshape(sigma_all, 1, 4);
sigma_p_shear(1, :) = reshape(sigma_p, 1, 4);
%sigma_w_shear(1, :) = reshape(sigma_w, 1, 4);
sigma_shear(1) = sigma_all(1, 2);
[~, Fval] = Energy_Disk_All_CL(xy_new, xc, yc, L0, N, D_all, Kp, Kpw, CL_pp, CL_pw);
%Fext_all(:, 1) = -Fval(2*N+1:2*N+8);

for it = 1:Nstep
    %fprintf('it: %d\n', it);
%     F_mat = [1, -it * dgamma;
%              0, 1];
    gamma_it = it * dgamma / 2;
    J = [1, gamma_it;
         gamma_it, 1];

    xy_c_strain = J * xy_c;
    xc_strain = xy_c_strain(1, :)';
    yc_strain = xy_c_strain(2, :)';
    [xy_new, ~] = FIRE_Disk_CL(xy_new, xc_strain, yc_strain, N, D_all, Fthresh, dt_fire, Nt_fire, Kp, Kpw, CL_pp, CL_pw);
    xy_t(:, it + 1) = [xy_new; xc_strain; yc_strain];
    gamma_all(it + 1) = it * dgamma;

    [sigma_all, sigma_p] = StressTensorAll_CL(xy_t(:, it + 1), L0, N, D_all, Kpw, CL_pw);
    sigma_all_shear(it + 1, :) = reshape(sigma_all, 1, 4);
    sigma_p_shear(it + 1, :) = reshape(sigma_p, 1, 4);
    %sigma_w_shear(it + 1, :) = reshape(sigma_w, 1, 4);
    sigma_shear(it + 1) = sigma_all(1, 2);
%     eigD = eig(sigma_tensor);
%     sigma_all(it) = 0.5 * abs(eigD(2) - eigD(1));
%     sigma_t = R_mat_stress * sigma_tensor * R_mat_stress';
%     sigma_all(it) = 0.5 * (sigma_t(1, 1) - sigma_t(2, 2));
%     sigma_tensor_all(it, :) = reshape(sigma_tensor, 1, 4);
    [~, Fval] = Energy_Disk_All_CL(xy_new, xc_strain, yc_strain, L0, N, D_all, Kp, Kpw, CL_pp, CL_pw);
%    Fext_all(:, it + 1) = -Fval(2*N+1:2*N+8);
end

G_fit_tot = polyfit(gamma_all, sigma_shear, 1);
G_fit_p = polyfit(gamma_all, sigma_p_shear(:, 2), 1);
%G_fit_w = polyfit(gamma_all, sigma_w_shear(:, 2), 1);
G_all = -[G_fit_tot(1); G_fit_p(1)];