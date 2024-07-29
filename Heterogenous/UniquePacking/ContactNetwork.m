function [ir, Z, CN, A] = ContactNetwork(x, y, D, N, gamma)
%%
ir = isRattler(x, y, D, N, gamma);
%%
R = D / 2;

L_all = [34.4; 50; 54; 50]; % wall lengths
[xc, yc] = GenerateCorners(L_all, gamma);

ift = [2; 3; 4; 1];
lx = xc(ift) - xc;
ly = yc(ift) - yc;
Cl_all = lx .* yc(ift) - ly .* xc(ift); % constant for line equation for four walls

Z_all = 0;
Z = zeros(N, 1);
CN = zeros(N * 6, 2);
A = zeros(N + 4);

for n = 1:N
    if ir(n)
        continue;
    end
    for w = 1:4
        % signed distance to Wall w, d > 0 if inside of the box
        if (lx(w) * y(n) - ly(w) * x(n) - Cl_all(w)) / L_all(w) < R 
            Z(n) = Z(n) + 1;
            Z_all = Z_all + 1;
            CN(Z_all, :) = [n, w + N];
            A(n, w + N) = 1;
        end
    end
    for m = n+1:N
        if ir(m)
            continue;
        end
        dx = x(m) - x(n);
        dy = y(m) - y(n);
        dr = sqrt(dx^2 + dy^2);
        if dr < D
            Z(n) = Z(n) + 1;
            Z(m) = Z(m) + 1;
            Z_all = Z_all + 1;
            CN(Z_all, :) = [n, m];
            A(n, m) = 1;
        end
    end
end

CN = CN(1:Z_all, :);
A = A + A';