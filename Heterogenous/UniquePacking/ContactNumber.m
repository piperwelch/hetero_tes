function [Z, Nc, Nnr] = ContactNumber(xy_all, D, N)
%%
xy = xy_all(1:2*N);
xc = xy_all(2*N+1:2*N+4);
yc = xy_all(2*N+5:2*N+8);

ift = [2; 3; 4; 1];
lx = xc(ift) - xc;
ly = yc(ift) - yc;
L_all = sqrt(lx.^2 + ly.^2);
Cl_all = lx .* yc(ift) - ly .* xc(ift); % constant for line equation for four walls
%%
x = xy(1:N);
y = xy(N+1:2*N);

ir = isRattler(x, y, xc, yc, D, N);
%%
Z = zeros(N, 1);
Nc = 0;

for n = 1:N
    if ir(n)
        continue;
    end
    for w = 1:4
        % signed distance to Wall w, d > 0 if inside of the box
        if (lx(w) * y(n) - ly(w) * x(n) - Cl_all(w)) / L_all(w) < D(n) / 2
            Z(n) = Z(n) + 1;
            Nc = Nc + 1;
        end
    end
    for m = n+1:N
        if ir(m)
            continue;
        end
        dx = x(m) - x(n);
        dy = y(m) - y(n);
        dr = sqrt(dx^2 + dy^2);
        if dr < 0.5 * (D(n) + D(m))
            Z(n) = Z(n) + 1;
            Z(m) = Z(m) + 1;
            Nc = Nc + 1;
        end
    end
end

Nnr = N - sum(double(ir));