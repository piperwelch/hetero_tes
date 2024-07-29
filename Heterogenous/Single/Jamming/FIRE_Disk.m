function [xy_equ, Eval] = FIRE_Disk(xy_all, N, D_all, L0, Kp, Kpw, Fthresh, dt_md, Nt)
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
a_fire = a_start;
delta_a_fire = 1 - a_fire;
dt = dt_md;
dt_half = dt / 2;
v_num = N*2;
Vel = zeros(v_num, 1);
[Eval, Acc] = Energy_Disk(xy_all, N, D_all, L0, Kp, Kpw);
if max(abs(Acc)) < Fthresh
    xy_equ = xy_all;
    return;
end

N_pp = 0; % number of P being positive
N_pn = 0; % number of P being negative
it = 0;
%% FIRE
while (true)
    it = it + 1;
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
            xy_all(1:N*2) = xy_all(1:N*2) - Vel * dt_half;
            Vel = zeros(v_num, 1);
        end
    end
    
    % MD using Verlet method
    Vel = Vel + Acc * dt_half;
    Vel = delta_a_fire * Vel + a_fire * norm(Vel) * Acc / norm(Acc);
    xy_all(1:N*2) = xy_all(1:N*2) + Vel * dt;
    [Eval, Acc] = Energy_Disk(xy_all, N, D_all, L0, Kp, Kpw);
    Vel = Vel + Acc * dt_half;
    
%     if mod(it, 10000) == 0
%         fprintf('it: %d  Max Force: %.5e\n', it, max(abs(Acc)));
%     end

    % exit when forces are smaller than threshold
    if max(abs(Acc)) < Fthresh || it > Nt
        break;
    end
end
%%
xy_equ = xy_all;