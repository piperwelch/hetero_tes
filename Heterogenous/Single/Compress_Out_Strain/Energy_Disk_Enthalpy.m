function [Eval, Fval] = Energy_Disk_Enthalpy(xy_p, xy_c_unit, N, D_all, L0, Kp, Kpw, Kw, VL)
%%
% Kp = 1; % particle-particle spring constant
% Kpw = 1; % particle-wall spring constant
% Kw = 1;
Np = N;
Np_all = 2 * Np;
Nc = 4;
Nc_all = 2 * Nc;
Ntot = Np_all + 1;

Eval = 0;
Fval = zeros(Np_all + 1, 1);
%% inter-particle energy
xp = xy_p(1:Np);
yp = xy_p(Np+1:Np_all);
for it = 1:size(VL, 1)
    i = VL(it, 1);
    j = VL(it, 2);
    dx = xp(j) - xp(i);
    Dnm = 0.5 * (D_all(i) + D_all(j));
    if abs(dx) < Dnm
        dy = yp(j) - yp(i);
        if abs(dy) < Dnm
            dr = sqrt(dx^2 + dy^2);
            if dr < Dnm
                dd = 1 - dr / Dnm;
                Eval = Eval + Kp * dd^2 / 2;
                
                F = Kp * dd / dr / Dnm;
                dFx = F * dx;
                dFy = F * dy;
                Fval(i) = Fval(i) - dFx;
                Fval(i + Np) = Fval(i + Np) - dFy;
                Fval(j) = Fval(j) + dFx;
                Fval(j + Np) = Fval(j + Np) + dFy;
            end
        end
    end
end
%% particle-wall energy
L = xy_p(Ntot);
xc = xy_c_unit(1:Nc) * L;
yc = xy_c_unit(Nc+1:Nc_all) * L;
ift = [2; 3; 4; 1];

A_box = 0.5 * sum(xc .* yc(ift) - xc(ift) .* yc);

for it = 1:4
    idx1 = it;
    idx2 = ift(it);
    lx = xc(idx2) - xc(idx1);
    ly = yc(idx2) - yc(idx1);
    lk_sq = lx^2 + ly^2;
    lk = sqrt(lk_sq);
    
    F = Kw * (L0 / lk - 1);
    Fval(Ntot) = Fval(Ntot) + F * lk_sq;
end

%% particle-wall energy
for i = 1:N
    R = D_all(i) / 2;
    x = xp(i);
    y = yp(i);
    for w = 1:4
        idx1 = w;
        idx2 = ift(w);
        lx = xc(idx2) - xc(idx1);
        ly = yc(idx2) - yc(idx1);
        Cl1 = lx * xc(idx1) + ly * yc(idx1);
        Cl2 = lx * xc(idx2) + ly * yc(idx2);
        if (lx * x + ly * y - Cl1) * (lx * x + ly * y - Cl2) > 0
            continue;
        end
        lk = sqrt(lx^2 + ly^2);
        Cl = lx * yc(idx2) - ly * xc(idx2);
        d = (lx * y - ly * x - Cl) / lk; % signed distance to Wall w, d > 0 if inside of the box
        if d < R
            dd = 1 - d / R;
            Eval = Eval + Kpw * dd^2 / 2;

            F = Kpw * dd / lk / R;
            Fx = F * ly;
            Fy = -F * lx;
            Fval(i) = Fval(i) - Fx;
            Fval(i + Np) = Fval(i + Np) - Fy;

            xcont = x + d * ly / lk;
            ycont = y - d * lx / lk;
            Fval(Ntot) = Fval(Ntot) + Fx * xcont + Fy * ycont;
        end
    end
end

Fval(Ntot) = Fval(Ntot) / A_box;