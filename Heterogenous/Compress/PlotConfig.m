function PlotConfig(xy_p_all, xy_c_all, N, Nvoxel_row, D_all, Wlist, linklist, ext_list, cn_on)
%%
Nvoxel = Nvoxel_row * Nvoxel_row;
Nvert = (Nvoxel_row + 1)^2;
Np = Nvoxel * N;
x = xy_p_all(1:Np);
y = xy_p_all(Np+1:2*Np);
xc = xy_c_all(1:Nvert);
yc = xy_c_all(Nvert+1:2*Nvert);

figure; hold on;

pcolor = [0 0.8 1]; % particle color
for k = 1:Np
    D = D_all(k);
    R = D / 2;
    rectangle('Position', [x(k) - R, y(k) - R, D, D], 'Curvature', [1 1],...
              'EdgeColor', pcolor, 'FaceColor', pcolor);
end

ift = [2; 3; 4; 1];
for i = 1:size(linklist, 1)
    j = linklist(i, 1);
    k = linklist(i, 2);
    plot([xc(j) xc(k)], [yc(j) yc(k)], 'k-', 'linewidth', 2);
end

for i = 1:Nvert
    plot(xc(i), yc(i), 'ro', 'linewidth', 1.5, 'markersize', 6, 'markerfacecolor', 'r')
end

for i = 1:size(ext_list, 1)
    wi = ext_list(i);
    plot(xc(wi), yc(wi), 'bo', 'linewidth', 1.5, 'markersize', 6, 'markerfacecolor', 'b')
end

x = reshape(x, N, Nvoxel);
y = reshape(y, N, Nvoxel);
D_all = reshape(D_all, N, Nvoxel);

if cn_on == 1 % plot contact network as well
    for nvoxel = 1:Nvoxel
        xv = x(:, nvoxel);
        yv = y(:, nvoxel);
        D_all_v = D_all(:, nvoxel);
        wlist = Wlist(:, nvoxel);
        xc_v = xc(wlist);
        yc_v = yc(wlist);

        for n = 1:N-1
            for m = n+1:N
                D = 0.5 * (D_all_v(n) + D_all_v(m));
                dx = xv(m) - xv(n);
                dy = yv(m) - yv(n);
                dr = sqrt(dx^2 + dy^2);
                if dr < D
                    plot([xv(n) xv(m)], [yv(n) yv(m)], 'r-', 'linewidth', 1.5);
                end
            end
        end

        for n = 1:N
            xp = xv(n);
            yp = yv(n);
            R = D_all_v(n) / 2;
            for w = 1:4
                idx1 = w;
                idx2 = ift(w);
                lx = xc_v(idx2) - xc_v(idx1);
                ly = yc_v(idx2) - yc_v(idx1);
                Cl1 = lx * xc_v(idx1) + ly * yc_v(idx1);
                Cl2 = lx * xc_v(idx2) + ly * yc_v(idx2);
                if (lx * xp + ly * yp - Cl1) * (lx * xp + ly * yp - Cl2) > 0
                    continue;
                end
                lk = sqrt(lx^2 + ly^2);
                Cl = lx * yc_v(idx2) - ly * xc_v(idx2);
                d = (lx * yp - ly * xp - Cl) / lk; % signed distance to Wall w, d > 0 if inside of the box
                if d < R
                    plot([xp xp + R * ly / lk], [yp yp - R * lx / lk], 'r-', 'linewidth', 1.5);
                end
            end
        end
    end
end

axis equal
axis off