function [Eval, Fval] = Energy_Disk(xy_all, N, Nvoxel, Nvert, D_all, L0, Wlist, linklist, Kp, Kpw, Kw)
%%
% Kp = 1; % particle-particle spring constant
% Kpw = 1; % particle-wall spring constant
% Kw = 1;
Np = N * Nvoxel;
Np_all = 2 * Np;
Nvert_all = 2 * Nvert;

Eval = 0;
Fp = zeros(Np_all, 1);
Fw = zeros(Nvert_all, 1);
%% inter-particle energy
xp = xy_all(1:Np);
yp = xy_all(Np+1:Np_all);
for nvoxel = 1:Nvoxel
    N_start = (nvoxel - 1) * N;
    for i = 1:N-1
        pi = i + N_start;
        for j = i+1:N
            pj = j + N_start;
            dx = xp(pj) - xp(pi);
            Dnm = 0.5 * (D_all(pi) + D_all(pj));
            if abs(dx) < Dnm
                
                dy = yp(pj) - yp(pi);
                if abs(dy) < Dnm
                    dr = sqrt(dx^2 + dy^2);
                    if dr < Dnm
                        dd = 1 - dr / Dnm;
                        Eval = Eval + Kp * dd^2 / 2;
                        
                        F = Kp * dd / dr / Dnm;
                        dFx = F * dx;
                        dFy = F * dy;
                        Fp(pi) = Fp(pi) - dFx;
                        Fp(pi + Np) = Fp(pi + Np) - dFy;
                        Fp(pj) = Fp(pj) + dFx;
                        Fp(pj + Np) = Fp(pj + Np) + dFy;
                    end
                end
            end
        end
    end
end
%% wall-wall energy
xc = xy_all(Np_all+1:Np_all+Nvert);
yc = xy_all(Np_all+Nvert+1:Np_all+Nvert_all);
ift = [2; 3; 4; 1];

for it = 1:size(linklist, 1)
    idx1 = linklist(it, 1);
    idx2 = linklist(it, 2);
    lx = xc(idx2) - xc(idx1);
    ly = yc(idx2) - yc(idx1);
    lk = sqrt(lx^2 + ly^2);
    
    Eval = Eval + Kw * (L0 - lk)^2 / 2;
    
    F = Kw * (L0 / lk - 1);
    Fx = F * lx;
    Fy = F * ly;
    Fw(idx1) = Fw(idx1) - Fx;
    Fw(idx1 + Nvert) = Fw(idx1 + Nvert) - Fy;
    Fw(idx2) = Fw(idx2) + Fx;
    Fw(idx2 + Nvert) = Fw(idx2 + Nvert) + Fy;
end

%% particle-wall energy
for nvoxel = 1:Nvoxel
    N_start = (nvoxel - 1) * N;
    wlist = Wlist(:, nvoxel);
    for i = 1:N
        pi = i + N_start;
        R = D_all(pi) / 2;
        x = xp(pi);
        y = yp(pi);
        for w = 1:4
            idx1 = wlist(w);
            idx2 = wlist(ift(w));
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
    %         fprintf('n: %d  w: %d  dnm: %.10e  dd: %.10e\n', i, w, d, R - d);
            if d < R
                dd = 1 - d / R;
                Eval = Eval + Kpw * dd^2 / 2;
    
                F = Kpw * dd / lk / R;
                Fp(pi) = Fp(pi) - F * ly;
                Fp(pi + Np) = Fp(pi + Np) + F * lx;
    
                dx = -d * lx / lk; % negative sign compared to DP code
                dy = -d * ly / lk;
                Fw(idx1) = Fw(idx1) + F * (yc(idx2) - y - dx);
                Fw(idx2) = Fw(idx2) - F * (yc(idx1) - y - dx);
                Fw(idx1 + Nvert) = Fw(idx1 + Nvert) - F * (xc(idx2) - x + dy);
                Fw(idx2 + Nvert) = Fw(idx2 + Nvert) + F * (xc(idx1) - x + dy);
            end
        end
    end
end
%%
Fval = [Fp; Fw];