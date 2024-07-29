function [xy_p_equ, xy_c_equ, Eval] = FIRE_Disk(xy_p_all, xy_c_all, N, Nvoxel, Nvert, D_all, L0,...
                                                Wlist, linklist, ext_list, Kp, Kpw, Kw, Fthresh, dt_md, Nt)
%% This is based on FIRE 2.0 provided in https://arxiv.org/pdf/1908.02038.pdf
%% Velocity verlet is chosen as the MD simulator
%% FIRE parameter for fast energy minimization
N_delay = 20;
N_pn_max = 2000;
f_inc = 1.1;
f_dec = 0.5;
a_start = 0.15;
f_a = 0.99;
% dt_md = 2e-3 / Epara.Kp;
dt_max = 10 * dt_md;
dt_min = 0.05 * dt_md;
% Nt = 2e5;
initialdelay = 1;
%% Initialization
Np = N * Nvoxel;
Np_tot = 2 * Np;
ext_list = [ext_list; ext_list + Nvert] + Np_tot;
xy_all = [xy_p_all; xy_c_all];
xy_c_fix = xy_all(ext_list);

a_fire = a_start;
delta_a_fire = 1 - a_fire;
dt = dt_md;
dt_half = dt / 2;
v_num = size(xy_all, 1);
Vel = zeros(v_num, 1);
[Eval, Acc] = Energy_Disk(xy_all, N, Nvoxel, Nvert, D_all, L0, Wlist, linklist, Kp, Kpw, Kw);
Acc(ext_list) = 0;
if max(abs(Acc)) < Fthresh
    xy_p_equ = xy_all(1:Np_tot);
    xy_c_equ = xy_all(Np_tot+1:Np_tot+2*Nvert);
    return;
end

N_pp = 0; % number of P being positive
N_pn = 0; % number of P being negative
%% FIRE
for it = 1:Nt
    % FIRE update
    P = Acc' * Vel;
    
    if P > 0
        N_pp = N_pp + 1;
        N_pn = 0;
        if N_pp > N_delay
            dt = min(f_inc * dt, dt_max);
            dt_half = dt / 2;
            a_fire = f_a * a_fire;
            delta_a_fire = 1 - a_fire;
        end
    else
        N_pp = 0;
        N_pn = N_pn + 1;
        if N_pn > N_pn_max
            break;
        end
        if (initialdelay < 0.5) || (it >= N_delay)
            if f_dec * dt > dt_min
                dt = f_dec * dt;
                dt_half = dt / 2;
            end
            a_fire = a_start;
            delta_a_fire = 1 - a_fire;
            xy_all = xy_all - Vel * dt_half;
            Vel = zeros(v_num, 1);
        end
    end
    
    % MD using Verlet method
    Vel = Vel + Acc * dt_half;
    Vel = delta_a_fire * Vel + a_fire * norm(Vel) * Acc / norm(Acc);
    xy_all = xy_all + Vel * dt;
    xy_all(ext_list) = xy_c_fix;
    [Eval, Acc] = Energy_Disk(xy_all, N, Nvoxel, Nvert, D_all, L0, Wlist, linklist, Kp, Kpw, Kw);
    Acc(ext_list) = 0;
    Vel = Vel + Acc * dt_half;
    
%     if mod(it, 10000) == 0
%         fprintf('it: %d  Max Force: %.5e\n', it, max(abs(Acc)));
%     end

    % exit when forces are smaller than threshold
    if max(abs(Acc)) < Fthresh
        break;
    end
end
%%
xy_p_equ = xy_all(1:Np_tot);
xy_c_equ = xy_all(Np_tot+1:Np_tot+2*Nvert);