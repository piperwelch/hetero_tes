function [xy_p_all_t, xy_c_all_t, L0_voxel_t, P_t] = Compress(xy_p_all, xy_c_all, N, Nvoxel_row, D_all, L0, L0_voxel, Kp, Kpw_ratio, Kw_ratio, test)
%% construct wall vertex list
Nvoxel = Nvoxel_row * Nvoxel_row;
Nvert = (Nvoxel_row + 1)^2;
Np = Nvoxel * N;
[Wlist, linklist, linklist_inner, ext_list] = ConstructTruss(Nvoxel_row);
%% FIRE parameters
rsc = 1.00000001;
Kw = Kw_ratio * Kp;
Kpw = Kpw_ratio * Kp;

dt_fire = 0.01 / sqrt(max([1, 4 * Kpw, Kw]));
Fthresh = 1E-14 * max([1, 4 * Kpw, Kw]);
Nt_fire = 1E7;

[~, P_old] = StressTensor([xy_p_all; xy_c_all], N, Nvoxel, Nvert, D_all, L0, Wlist, linklist, Kp, Kpw, Kw);

xy_p_all_t = zeros(2 * Np, 1000);
xy_c_all_t = zeros(2 * Nvert, 1000);
L0_voxel_t = zeros(1000, 1);
P_t = zeros(1000, 1);
count = 1;
xy_p_all_t(:, 1) = xy_p_all;
xy_c_all_t(:, 1) = xy_c_all;
L0_voxel_t(1) = L0_voxel;
P_t(1) = P_old;

L0_voxel = L0_voxel / rsc;
xy_p_all = xy_p_all / rsc;
xy_c_all = xy_c_all / rsc;
%%
while true
    [xy_p_all, xy_c_all, ~] = FIRE_Disk(xy_p_all, xy_c_all, N, Nvoxel, Nvert, D_all, L0,...
                                        Wlist, linklist, ext_list, Kp, Kpw, Kw, Fthresh, dt_fire, Nt_fire);
    [~, P] = StressTensor([xy_p_all; xy_c_all], N, Nvoxel, Nvert, D_all, L0, Wlist, linklist, Kp, Kpw, Kw);
    if test == 1
        fprintf("P: %.5e  L: %.5e\n", P, L0_voxel);
    end

    if P < P_old || P < 1e-7
        L0_voxel = L0_voxel / rsc;
        xy_p_all = xy_p_all / rsc;
        xy_c_all = xy_c_all / rsc;
        continue;
    end

    if P > 1e-2
        break;
    end

    if P < 1e-6 && P >= 1e-7
        rsc = 1.00000001;
    elseif P < 1e-5 && P >= 1e-6
        rsc = 1.0000001;
    elseif P < 1e-4 && P >= 1e-5
        rsc = 1.000001;
    elseif P < 1e-3 && P >= 1e-4
        rsc = 1.00001;
    elseif P < 1e-2 && P >= 1e-3
        rsc = 1.0001;
    end

    P_old = P;
    count = count + 1;
    xy_p_all_t(:, count) = xy_p_all;
    xy_c_all_t(:, count) = xy_c_all;
    L0_voxel_t(count) = L0_voxel;
    P_t(count) = P;
    L0_voxel = L0_voxel / rsc;
    xy_p_all = xy_p_all / rsc;
    xy_c_all = xy_c_all / rsc;
end
%%
xy_p_all_t = xy_p_all_t(:, 1:count);
xy_c_all_t = xy_c_all_t(:, 1:count);
L0_voxel_t = L0_voxel_t(1:count);
P_t = P_t(1:count);