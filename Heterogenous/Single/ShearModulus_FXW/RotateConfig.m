function xy_new = RotateConfig(N, xy_all, ridx1, theta)
%%
x = xy_all(1:N);
y = xy_all(N+1:2*N);
xc = xy_all(2*N+1:2*N+4);
yc = xy_all(2*N+5:2*N+8);

R_mat = [cos(theta), -sin(theta);
         sin(theta), cos(theta)];

xy_rot = R_mat * [x' - xc(ridx1), xc' - xc(ridx1); y' - yc(ridx1), yc' - yc(ridx1)];
x = xy_rot(1, 1:N)';
y = xy_rot(2, 1:N)';
xc = xy_rot(1, N+1:N+4)';
yc = xy_rot(2, N+1:N+4)';

xy_new = [x; y; xc; yc];