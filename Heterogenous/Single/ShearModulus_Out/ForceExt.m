function Fw = ForceExt(xy_p, xc, yc, L_all, N, D_all, Fiso)
%%
% Kp = 1; % particle-particle spring constant
Kpw = 1; % particle-wall spring constant
Kw = 1;
Np = N;
Np_all = 2 * Np;
Nc = 4;
Nc_all = 2 * Nc;

xp = xy_p(1:Np);
yp = xy_p(Np+1:Np_all);

Fw = zeros(Nc_all, 1);
%% wall-wall energy
ift = [2; 3; 4; 1];

for it = 1:4
    idx1 = it;
    idx2 = ift(it);
    lx = xc(idx2) - xc(idx1);
    ly = yc(idx2) - yc(idx1);
    lk = sqrt(lx^2 + ly^2);
    
    F = Kw * (L_all(it) / lk - 1);
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
            F = Kpw * dd / lk;

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
jft = [4; 1; 2; 3];
lx = xc(ift) - xc;
ly = yc(ift) - yc;
Fw = Fiso * [-ly - ly(jft), lx + lx(jft)] / 2 - [Fw(1:Nc), Fw(Nc+1:Nc_all)];