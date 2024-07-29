function [sigma_all, sigma_p, sigma_w] = StressTensorAll(xy, L0, N, D_all, Kw)
%%
Kp = 1; % particle-particle spring constant
Kpw = 1; % particle-wall spring constant
% Kw = 1;
Np = N;
Np_all = 2 * Np;
Nc = 4;
Nc_all = 2 * Nc;

sigma_p = zeros(2);
sigma_w = zeros(2);
%% inter-particle energy
xp = xy(1:Np);
yp = xy(Np+1:Np_all);
% for i = 1:N-1
%     for j = i+1:N
%         dx = xp(j) - xp(i);
%         Dnm = 0.5 * (D_all(i) + D_all(j));
%         if abs(dx) < Dnm
%             dy = yp(j) - yp(i);
%             if abs(dy) < Dnm
%                 dr = sqrt(dx^2 + dy^2);
%                 if dr < Dnm
%                     dd = Dnm - dr;
%                     Eval = Eval + Kp * dd^2 / 2;
%                     
%                     F = Kp * dd / dr;
%                     dFx = F * dx;
%                     dFy = F * dy;
%                     Fp(i) = Fp(i) - dFx;
%                     Fp(i + Np) = Fp(i + Np) - dFy;
%                     Fp(j) = Fp(j) + dFx;
%                     Fp(j + Np) = Fp(j + Np) + dFy;
%                 end
%             end
%         end
%     end
% end
%% wall-wall energy
xc = xy(Np_all+1:Np_all+Nc);
yc = xy(Np_all+Nc+1:Np_all+Nc_all);
ift = [2; 3; 4; 1];

A_box = 0.5 * sum(xc .* yc(ift) - xc(ift) .* yc);

for it = 1:4
    idx1 = it;
    idx2 = ift(it);
    lx = xc(idx2) - xc(idx1);
    ly = yc(idx2) - yc(idx1);
    lk = sqrt(lx^2 + ly^2);
    
    F = Kw * (L0 / lk - 1);
    Fx = F * lx;
    Fy = F * ly;
    sigma_w(1, 1) = sigma_w(1, 1) + Fx * lx;
    sigma_w(1, 2) = sigma_w(1, 2) + Fx * ly;
    sigma_w(2, 1) = sigma_w(2, 1) + Fy * lx;
    sigma_w(2, 2) = sigma_w(2, 2) + Fy * ly;
end

sigma_w = sigma_w / A_box;
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
            F = Kpw * dd / lk;
            Fx = F * ly;
            Fy = -F * lx;
            xcont = x + d * ly / lk;
            ycont = y - d * lx / lk;

            sigma_p(1, 1) = sigma_p(1, 1) + Fx * xcont;
            sigma_p(1, 2) = sigma_p(1, 2) + Fx * ycont;
            sigma_p(2, 1) = sigma_p(2, 1) + Fy * xcont;
            sigma_p(2, 2) = sigma_p(2, 2) + Fy * ycont;
        end
    end
end
sigma_p = sigma_p / A_box;
%%
sigma_all = sigma_w + sigma_p;