function [xy_p_new, xy_c_new] = RotateConfig(xy_p_all, xy_c_all, Np, Nvert, theta)
%%
xp = xy_p_all(1:Np);
yp = xy_p_all(Np+1:2*Np);
xc = xy_c_all(1:Nvert);
yc = xy_c_all(Nvert+1:2*Nvert);

R_mat = [cos(theta), -sin(theta);
         sin(theta), cos(theta)];

xy_rot = R_mat * [xp', xc'; yp', yc'];
xp = xy_rot(1, 1:Np)';
yp = xy_rot(2, 1:Np)';
xc = xy_rot(1, Np+1:Np+Nvert)';
yc = xy_rot(2, Np+1:Np+Nvert)';

xy_p_new = [xp; yp]; 
xy_c_new = [xc; yc];