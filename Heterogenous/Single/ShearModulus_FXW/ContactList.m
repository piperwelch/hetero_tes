function [CL_pp, CL_pw] = ContactList(xy, N, D_all)
%%
Np = N;
Np_all = 2 * Np;
Nc = 4;
Nc_all = 2 * Nc;

CL_pp = zeros(6 * N, 2);
CL_pw = zeros(6 * N, 2);
count_pp = 0;
count_pw = 0;
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
                    count_pp = count_pp + 1;
                    CL_pp(count_pp, :) = [i, j];
                end
            end
        end
    end
end
%% wall-wall energy
xc = xy(Np_all+1:Np_all+Nc);
yc = xy(Np_all+Nc+1:Np_all+Nc_all);
ift = [2; 3; 4; 1];

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
            count_pw = count_pw + 1;
            CL_pw(count_pw, :) = [i, w];
        end
    end
end
%%
CL_pp = CL_pp(1:count_pp, :);
CL_pw = CL_pw(1:count_pw, :);