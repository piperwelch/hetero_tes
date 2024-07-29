function [sigma, P] = StressTensor(xy_all, N, D_all, Kp, Kpw)
%%
xy = xy_all(1:2*N);
xy_c = xy_all(2*N+1:2*N+8);
sigma = zeros(2);
%% inter-particle energy
for n = 1:N-1
    for m = n+1:N
        Dnm = 0.5 * (D_all(n) + D_all(m));
        dx = xy(m) - xy(n);
        dy = xy(m + N) - xy(n + N);
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

%% particle-wall energy
ift = [2; 3; 4; 1];
xc = xy_c(1:4);
yc = xy_c(5:8);
lx = xc(ift) - xc;
ly = yc(ift) - yc;
L_all = sqrt(lx.^2 + ly.^2);
Cl_all = lx .* yc(ift) - ly .* xc(ift);

for n = 1:N
    R = D_all(n) / 2;
    ny = n + N;
    x = xy(n);
    y = xy(ny);
    for w = 1:4
        d = (lx(w) * y - ly(w) * x - Cl_all(w)) / L_all(w); % signed distance to Wall w, d > 0 if inside of the box
        if d < R    
            F = Kpw * (R - d) / L_all(w) / R^2;
            dFx = F * ly(w);
            dFy = -F * lx(w);
            dx = d * ly(w) / L_all(w);
            dy = -d * lx(w) / L_all(w);
            sigma(1, 1) = sigma(1, 1) + dFx * dx;
            sigma(1, 2) = sigma(1, 2) + dFx * dy;
            sigma(2, 1) = sigma(2, 1) + dFy * dx;
            sigma(2, 2) = sigma(2, 2) + dFy * dy;
        end
    end
end
%%
A_box = 0.5 * sum(xc .* yc(ift) - xc(ift) .* yc);
sigma = sigma / A_box;
P = 0.5 * (sigma(1, 1) + sigma(2, 2));