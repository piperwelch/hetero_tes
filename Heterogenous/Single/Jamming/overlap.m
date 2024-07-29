function ol = overlap(xy, N, D_all, xy_c)
%%
ol = 0;

for n = 1:N-1
    for m = n+1:N
        Dnm = 0.5 * (D_all(n) + D_all(m));
        dx = xy(m) - xy(n);
        dy = xy(m + N) - xy(n + N);
        if sqrt(dx^2 + dy^2) < Dnm
            ol = 1;
            return;
        end
    end
end

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
            ol = 1;
            return;
        end
    end
end