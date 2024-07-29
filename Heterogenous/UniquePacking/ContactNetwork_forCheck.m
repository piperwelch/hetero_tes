function [Nnr, Nnr_b, Nnr_s, Z_all, A] = ContactNetwork_forCheck(x, y, xc, yc, D, N)
%%
D_mean = mean(D);
ir = isRattler(x, y, xc, yc, D, N);
x(ir) = [];
y(ir) = [];
D(ir) = [];
Nnr = size(x, 1);
Nnr_b = size(D(D > D_mean), 1);
Nnr_s = size(D(D < D_mean), 1);
%%
ift = [2; 3; 4; 1];
lx = xc(ift) - xc;
ly = yc(ift) - yc;
L_all = sqrt(lx.^2 + ly.^2);
Cl_all = lx .* yc(ift) - ly .* xc(ift); % constant for line equation for four walls

Z_all = 0;
% Z = zeros(Nnr, 1);
% CN = zeros(Nnr * 6, 2);
A = zeros(Nnr + 4);

for n = 1:Nnr
    for w = 1:4
        % signed distance to Wall w, d > 0 if inside of the box
        if (lx(w) * y(n) - ly(w) * x(n) - Cl_all(w)) / L_all(w) < D(n) / 2
%             Z(n) = Z(n) + 1;
            Z_all = Z_all + 1;
%             CN(Z_all, :) = [n, w + Nnr];
            A(n, w + Nnr) = 1;
        end
    end
    for m = n+1:Nnr
        dx = x(m) - x(n);
        dy = y(m) - y(n);
        dr = sqrt(dx^2 + dy^2);
        if dr < 0.5 * (D(n) + D(m))
%             Z(n) = Z(n) + 1;
%             Z(m) = Z(m) + 1;
            Z_all = Z_all + 1;
%             CN(Z_all, :) = [n, m];
            A(n, m) = 1;
        end
    end
end

% CN = CN(1:Z_all, :);
A = A + A';