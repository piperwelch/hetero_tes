function [xy_p_equ, xy_c_equ, Eval] = FIRE_Disk_Enthalpy(xy_p, xy_c, N, D_all, L0, Fthresh, dt_md, Nt, Kp, Kpw, Kw, Pt)
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
a_fire = a_start;
delta_a_fire = 1 - a_fire;
dt = dt_md;
dt_half = dt / 2;
%% VL initialization
xy_save = xy_p;
VL = zeros(N * 10, 2);
[VL, xy_save] = VerletList_Disk(xy_p, xy_save, N, D_all, 1, VL);
%% Initialization
xc = xy_c(1:4);
yc = xy_c(5:8);
ift = [2; 3; 4; 1];
L = sqrt((xc(2) - xc(1))^2 + (yc(2) - yc(1))^2);
A_coef = 0.5 * sum(xc .* yc(ift) - xc(ift) .* yc) / L^2;
xy_c_unit = xy_c / L;
xy_p = [xy_p; L];
m_wall = sqrt(N); % for stable energy minimization

v_num = size(xy_p, 1);
Vel = zeros(v_num, 1);
[Eval, Acc] = Energy_Disk_Enthalpy(xy_p, xy_c_unit, N, D_all, L0, Kp, Kpw, Kw, VL);
Acc(v_num) = (Acc(v_num) - 2 * Pt) * xy_p(v_num) * A_coef / m_wall;

if max(abs(Acc(1:v_num-1))) < Fthresh && abs(Acc(v_num)) < 2e-5 * A_coef * Pt * xy_p(v_num)
    xy_p_equ = xy_p;
    xy_c_equ = xy_c;
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
            xy_p = xy_p - Vel * dt_half;
            Vel = zeros(v_num, 1);
        end
    end
    
    % MD using Verlet method
    Vel = Vel + Acc * dt_half;
    Vel = delta_a_fire * Vel + a_fire * norm(Vel) * Acc / norm(Acc);
    xy_p = xy_p + Vel * dt;
    [VL, xy_save] = VerletList_Disk(xy_p, xy_save, N, D_all, 0, VL);
    [Eval, Acc] = Energy_Disk_Enthalpy(xy_p, xy_c_unit, N, D_all, L0, Kp, Kpw, Kw, VL);
    Acc(v_num) = (Acc(v_num) - 2 * Pt) * xy_p(v_num) * A_coef / m_wall;
    Vel = Vel + Acc * dt_half;
    
    % For recording
%     Eval(it) = E;
%     F2norm(it) = sum(grad.^2);

    % exit when forces are smaller than threshold
    if max(abs(Acc(1:v_num-1))) < Fthresh && abs(Acc(v_num)) < 2e-5 * A_coef * Pt * xy_p(v_num)
        break;
    end
end
%%
xy_p_equ = xy_p(1:2*N);
xy_c_equ = xy_c_unit * xy_p(end);