function [sigma, P] = StressTensor(xy_all, N, Nvoxel, Nvert, D_all, L0, Wlist, linklist, Kp, Kpw, Kw)
%%
% Kp = 1;
% Kpw = 1;
% Kw = 1;
Np = N * Nvoxel;
Np_all = 2 * Np;
Nvert_all = 2 * Nvert;
xp = xy_all(1:Np);
yp = xy_all(Np+1:Np_all);
sigma = zeros(2);
%% inter-particle energy
for nvoxel = 1:Nvoxel
    N_start = (nvoxel - 1) * N;
    for n = 1:N-1
        pn = n + N_start;
        for m = n+1:N
            pm = m + N_start;
            Dnm = 0.5 * (D_all(pn) + D_all(pm));
            dx = xp(pm) - xp(pn);
            dy = yp(pm) - yp(pn);
            dr = sqrt(dx^2 + dy^2);
            if dr < Dnm
                F = Kp * (Dnm / dr - 1) / Dnm^2;
                dFx = F * dx;
                dFy = F * dy;
                sigma(1, 1) = sigma(1, 1) + dFx * dx;
                sigma(1, 2) = sigma(1, 2) + dFx * dy;
                sigma(2, 1) = sigma(2, 1) + dFy * dx;
                sigma(2, 2) = sigma(2, 2) + dFy * dy;
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
    
    F = Kw * (L0 / lk - 1);
    dFx = F * lx;
    dFy = F * ly;
    sigma(1, 1) = sigma(1, 1) + dFx * lx;
    sigma(1, 2) = sigma(1, 2) + dFx * ly;
    sigma(2, 1) = sigma(2, 1) + dFy * lx;
    sigma(2, 2) = sigma(2, 2) + dFy * ly;
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
                F = Kpw * (R - d) / lk / R^2;
                dFx = F * ly;
                dFy = -F * lx;
                dx = d * ly / lk;
                dy = -d * lx / lk;
                sigma(1, 1) = sigma(1, 1) + dFx * dx;
                sigma(1, 2) = sigma(1, 2) + dFx * dy;
                sigma(2, 1) = sigma(2, 1) + dFy * dx;
                sigma(2, 2) = sigma(2, 2) + dFy * dy;
            end
        end
    end
end

%%
Nvoxel_row = sqrt(Nvoxel);
La = [xc(Nvoxel_row + 1) - xc(1); yc(Nvoxel_row + 1) - yc(1)];
Lb = [xc(Nvoxel_row * (Nvoxel_row + 1) + 1) - xc(1); yc(Nvoxel_row * (Nvoxel_row + 1) + 1) - yc(1)];
A_box = abs(La(1) * Lb(2) - La(2) * Lb(1));
sigma = sigma / A_box;
P = 0.5 * (sigma(1, 1) + sigma(2, 2));