function xy_new = RotateConfig_Voxel(N, xy_all, ridx1, ridx2, vertex_reorder)
%%
x = xy_all(1:N);
y = xy_all(N+1:2*N);
xc = xy_all(2*N+1:2*N+4);
yc = xy_all(2*N+5:2*N+8);

rt = atan2(yc(ridx2) - yc(ridx1), xc(ridx2) - xc(ridx1));
R_mat = [cos(rt), sin(rt); -sin(rt), cos(rt)];

xy_rot = R_mat * [x - xc(ridx1), xc - xc(ridx1); y - yc(ridx1), yc - yc(ridx1)];
x = xy_rot(1, 1:N)';
y = xy_rot(2, 1:N)';
xc = xy_rot(1, N+1:N+4)';
yc = xy_rot(2, N+1:N+4)';
xc = xc(vertex_reorder);
yc = yc(vertex_reorder);

xy_new = [x; y; xc; yc];