function [xy_new, theta] = RotateConfig_Initial(N, xy_all)
%%
x = xy_all(1:N);
y = xy_all(N+1:2*N);
xc = xy_all(2*N+1:2*N+4);
yc = xy_all(2*N+5:2*N+8);

dr1 = [xc(2) - xc(1), yc(2) - yc(1)];
dr2 = [xc(4) - xc(1), yc(4) - yc(1)];
theta = acos(dr1 * dr2' / norm(dr1) / norm(dr2));

ridx1 = 1;
ridx2 = 2;
rt = atan2(yc(ridx2) - yc(ridx1), xc(ridx2) - xc(ridx1));
R_mat = [cos(rt), sin(rt); -sin(rt), cos(rt)];

xy_rot = R_mat * [x - xc(ridx1), xc - xc(ridx1); y - yc(ridx1), yc - yc(ridx1)];
x = xy_rot(1, 1:N)';
y = xy_rot(2, 1:N)';
xc = xy_rot(1, N+1:N+4)';
yc = xy_rot(2, N+1:N+4)';

xy_new = [x; y; xc; yc];