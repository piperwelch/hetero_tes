function [xy_all_t, P_t] = Compress_Enthalpy(xy_all, L0, N, D_all, Kpw_ratio, Kw_ratio, test)
%%
xy_all_t = zeros(2 * N + 8, 101);
P_t = logspace(-7, -2, 101)';

Kp = 1;
Kpw = Kpw_ratio * Kp;
Kw = Kw_ratio * Kp;
Fthresh = 1E-14 * max([1, 4 * Kpw, Kw]);
dt_fire = 0.01 / max([1, 4 * Kpw, Kw]);
Nt_fire = 1E7;

Np_tot = 2 * N;
xy_all = RotateConfig_Initial(N, xy_all, 1, 2);
xy_p = xy_all(1:Np_tot);
xy_c = xy_all(Np_tot+1:Np_tot+8);
%%
for it = 1:101
    if test == 1
        fprintf('it: %d\n', it);
    end
    [xy_p, xy_c, ~] = FIRE_Disk_Enthalpy(xy_p, xy_c, N, D_all, L0, Fthresh, dt_fire, Nt_fire, Kp, Kpw, Kw, P_t(it)); % energy minimization, equivalently force balance
    xy_all_t(:, it) = [xy_p; xy_c];
end