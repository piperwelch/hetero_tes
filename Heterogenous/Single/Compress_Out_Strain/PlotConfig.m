function PlotConfig(xy_all, N, D_all, cn_on)
%%
Np = N;
Np_tot = 2 * Np;
Nc = 4;
Nc_tot = 2 * Nc;
x = xy_all(1:Np);
y = xy_all(Np+1:Np_tot);
xc = xy_all(Np_tot+1:Np_tot+Nc);
yc = xy_all(Np_tot+Nc+1:Np_tot+Nc_tot);

figure; hold on;

pcolor = [0 0.8 1]; % particle color
for k = 1:N
    D = D_all(k);
    R = D / 2;
    rectangle('Position', [x(k) - R, y(k) - R, D, D], 'Curvature', [1 1],...
              'EdgeColor', pcolor, 'FaceColor', pcolor);
end

ift = [2; 3; 4; 1];
for i = 1:4
    plot([xc(i) xc(ift(i))], [yc(i) yc(ift(i))],...
         'ko-', 'linewidth', 1.5, 'markersize', 4, 'markerfacecolor', 'k');
end

if cn_on == 1 % plot contact network as well
    for n = 1:N-1
        for m = n+1:N
            D = 0.5 * (D_all(n) + D_all(m));
            dx = x(m) - x(n);
            dy = y(m) - y(n);
            dr = sqrt(dx^2 + dy^2);
            if dr < D
                plot([x(n) x(m)], [y(n) y(m)], 'r-', 'linewidth', 1.5);
            end
        end
    end

    for i = 1:N
        xp = x(i);
        yp = y(i);
        R = D_all(i) / 2;
        for w = 1:4
            idx1 = w;
            idx2 = ift(w);
            lx = xc(idx2) - xc(idx1);
            ly = yc(idx2) - yc(idx1);
            Cl1 = lx * xc(idx1) + ly * yc(idx1);
            Cl2 = lx * xc(idx2) + ly * yc(idx2);
            if (lx * xp + ly * yp - Cl1) * (lx * xp + ly * yp - Cl2) > 0
                continue;
            end
            lk = sqrt(lx^2 + ly^2);
            Cl = lx * yc(idx2) - ly * xc(idx2);
            d = (lx * yp - ly * xp - Cl) / lk; % signed distance to Wall w, d > 0 if inside of the box
            if d < R
                plot([xp xp + R * ly / lk], [yp yp - R * lx / lk], 'r-', 'linewidth', 1.5);
            end
        end
    end
end

axis equal
axis off