function [Eval, Fval] = Energy_Disk(xy, L0, N, D_all, Kw)
%%
Kp = 1; % particle-particle spring constant
Kpw = 1; % particle-wall spring constant
% Kw = 1;
Np = N;
Np_all = 2 * Np;
Nc = 4;
Nc_all = 2 * Nc;

Eval = 0;
Fp = zeros(Np_all, 1);
Fw = zeros(Nc_all, 1);
%% inter-particle energy
xp = xy(1:Np);
yp = xy(Np+1:Np_all);
for i = 1:N-1
    for j = i+1:N
        dx = xp(j) - xp(i);
        Dnm = 0.5 * (D_all(i) + D_all(j));
        if abs(dx) < Dnm
            dy = yp(j) - yp(i);
            if abs(dy) < Dnm
                dr = sqrt(dx^2 + dy^2);
                if dr < Dnm
                    dd = Dnm - dr;
                    Eval = Eval + Kp * dd^2 / 2;
                    
                    F = Kp * dd / dr;
                    dFx = F * dx;
                    dFy = F * dy;
                    Fp(i) = Fp(i) - dFx;
                    Fp(i + Np) = Fp(i + Np) - dFy;
                    Fp(j) = Fp(j) + dFx;
                    Fp(j + Np) = Fp(j + Np) + dFy;
                end
            end
        end
    end
end
%% wall-wall energy
xc = xy(Np_all+1:Np_all+Nc);
yc = xy(Np_all+Nc+1:Np_all+Nc_all);
ift = [2; 3; 4; 1];

for it = 1:4
    idx1 = it;
    idx2 = ift(it);
    lx = xc(idx2) - xc(idx1);
    ly = yc(idx2) - yc(idx1);
    lk = sqrt(lx^2 + ly^2);
    
    Eval = Eval + Kw * (L0 - lk)^2 / 2;
    
    F = Kw * (L0 / lk - 1);
    Fx = F * lx;
    Fy = F * ly;
    Fw(idx1) = Fw(idx1) - Fx;
    Fw(idx1 + Nc) = Fw(idx1 + Nc) - Fy;
    Fw(idx2) = Fw(idx2) + Fx;
    Fw(idx2 + Nc) = Fw(idx2 + Nc) + Fy;
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
            dd = R - d;
            Eval = Eval + Kpw * dd^2 / 2;
            F = Kpw * dd / lk;
            Fp(i) = Fp(i) - F * ly;
            Fp(i + Np) = Fp(i + Np) + F * lx;

            dx = -d * lx / lk; % negative sign compared to DP code
            dy = -d * ly / lk;
            Fw(idx1) = Fw(idx1) + F * (yc(idx2) - y - dx);
            Fw(idx2) = Fw(idx2) - F * (yc(idx1) - y - dx);
            Fw(idx1 + Nc) = Fw(idx1 + Nc) - F * (xc(idx2) - x + dy);
            Fw(idx2 + Nc) = Fw(idx2 + Nc) + F * (xc(idx1) - x + dy);
        end
    end
end
%%
Fval = [Fp; Fw];