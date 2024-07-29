function [xy_p_all, xy_c_all, L0_voxel, D_all_new] = Jamming(Ncons, xy_voxel, N, D_all, L0, Nvoxel_row, Kp, Kw_ratio, Kpw_ratio, test)

%{ Generates jammed configs 

%}

%% construct wall vertex list
Nvoxel = Nvoxel_row * Nvoxel_row;
Nvert = (Nvoxel_row + 1)^2;
Np = Nvoxel * N;
[Wlist, linklist, ~, ext_list] = ConstructTruss(Nvoxel_row);

%% replicate particle positions
theta = ones(Ncons, 1);
xy_voxel_new = ones(Ncons, 2*N + 8);

for i = 1:Ncons
    [xy_voxel_new(i, :), theta(i)] = RotateConfig_Initial(N, xy_voxel(i,:));
end

xy_p_new = xy_voxel_new(:, 1:2*N);
La = [L0; 0];
Lb = [cos(theta(1)); sin(theta(1))] * L0;

D_all_new = zeros(Nvoxel * 4, 1);

xy_p_all = zeros(2 * Np, 1);
count = 0;
for i = 1:Nvoxel_row
    for j = 1:Nvoxel_row
        N_start = count * N;
        count = count + 1;
        conid = 1 + rem(i + j, 2);
        D_all_new((count-1) * 4 + 1 : count*4) = D_all(conid, :);
        for n = 1:N         
            xy_p_all(n + N_start) = xy_p_new(conid, n) + (j - 1) * La(1) + (i - 1) * Lb(1);
            xy_p_all(n + N_start + Np) = xy_p_new(conid, n + N) + (j - 1) * La(2) + (i - 1) * Lb(2);
        end
    end
end

xy_c_all = zeros(2 * Nvert, 1);
count = 0;
for i = 1:Nvoxel_row+1
    for j = 1:Nvoxel_row+1
        count = count + 1;
        xy_c_all(count) = (j - 1) * La(1) + (i - 1) * Lb(1);
        xy_c_all(count + Nvert) = (j - 1) * La(2) + (i - 1) * Lb(2);
    end
end

%% FIRE parameters
rsc = 1.000001;
Pt_l = 1E-7; % low bound of total energy for the system to be deemed as jammed
Pt_h = 1.01 * Pt_l; % high bound of total energy for the system to be deemed as jammed

Kpw = Kpw_ratio * Kp;
Kw = Kw_ratio * Kp;

dt_fire = 0.01 / sqrt(max([1, 4 * Kpw, Kw]));
Fthresh = 1E-14 * max([1, 4 * Kpw, Kw]);
Nt_fire = 1E7;

[xy_p_all, xy_c_all, ~] = FIRE_Disk(xy_p_all, xy_c_all, N, Nvoxel, Nvert, D_all_new, L0,...
                                    Wlist, linklist, ext_list, Kp, Kpw, Kw, Fthresh, dt_fire, Nt_fire);
[~, P] = StressTensor([xy_p_all; xy_c_all], N, Nvoxel, Nvert, D_all_new, L0, Wlist, linklist, Kp, Kpw, Kw);

L0_voxel = L0;
if P > Pt_h
    L_l = L0_voxel;
    L_h = -1.0;
    L0_voxel = L0_voxel * 1.00001;
    xy_p_all = xy_p_all * 1.00001;
    xy_c_all = xy_c_all * 1.00001;
elseif P < Pt_l
    xy_p_old = xy_p_all;
    xy_c_old = xy_c_all;
    L_h = L0_voxel;
    L_l = -1.0;
    L0_voxel = L0_voxel / 1.00001;
    xy_p_all = xy_p_all / 1.00001;
    xy_c_all = xy_c_all / 1.00001;
else
    return;
end
%% 
while true
    [xy_p_all, xy_c_all, ~] = FIRE_Disk(xy_p_all, xy_c_all, N, Nvoxel, Nvert, D_all_new, L0,...
                                    Wlist, linklist, ext_list, Kp, Kpw, Kw, Fthresh, dt_fire, Nt_fire);
    [~, P] = StressTensor([xy_p_all; xy_c_all], N, Nvoxel, Nvert, D_all_new, L0, Wlist, linklist, Kp, Kpw, Kw);
    %PlotConfig(xy_p_all, xy_c_all, N, Nvoxel_row, D_all_new, Wlist, linklist, ext_list, 1)
    
    if test == 1
        fprintf("P: %.5e   L: %.5e\n", P, L0_voxel);
    end

    if P < Pt_l
        xy_p_old = xy_p_all;
        xy_c_old = xy_c_all;
        L_h = L0_voxel;
        if L_l > 0.0
            L0_voxel = (L_h + L_l) / 2.0;
            L_l = -1.0;
        else
            L0_voxel = L0_voxel / rsc;
        end
        xy_p_all = xy_p_old * L0_voxel / L_h;
        xy_c_all = xy_c_old * L0_voxel / L_h;
    elseif P > Pt_h
        L_l = L0_voxel;
        L0_voxel = (L_h + L_l) / 2.0;
        xy_p_all = xy_p_old * L0_voxel / L_h;
        xy_c_all = xy_c_old * L0_voxel / L_h;
    else
        break;
    end
    if ((L_l > 0.0) && (abs(L_h / L_l - 1.0) < 1E-14))
        break;
    end
end