function [xy_all, D_all, D] = Jamming(N, pid, ratio, Dtype, Kp, Kw_ratio, Kpw_ratio, test, theta)

%{
Generates jammed configs in boundaries at fixed theta (rhombus with length L)

input:
------
N (sc): particle number
pid (int): random seed for particle initialization 
ratio (sc): bidispersity ratio for particle sizes (radius big / radius small)
Dtype (0 or 1): for odd N, if 1 then more small particles; if 0 then more big particles
Kp (sc): particle-particle stiffness
Kw_ratio (sc): ratio between wall and p-p stiffness
Kpw_ratio (sc): ratio betwteen p-w and p-p stiffness
test (0 or 1): if 1, print intermediate values
theta (sc): angle in radians for box shape
-----

output:
------
xy_all [1 x (2N + 8)]: particle positions and wall positions, arranged like [x1 x2 x3.. xN y1 .. yN || vx1 ... vx4 vy1 ... vy4] where particles are first, and then box positions 
D_all (Nx1): actual particle diameters
D (Nx1): smaller particles only
------
%}

%%
L0 = 1;
Nvert = 4; %always true for single voxel
Nh = floor(N / 2);

%set two sizes of particles
if mod(N, 2) ~= 0
    if Dtype == 1 % 1 more small particle than big ones
        Gn = [ratio * ones(Nh, 1); ones(N - Nh, 1)];
    else % 1 more big particle than small ones
        Gn = [ratio * ones(Nh + 1, 1); ones(N - Nh - 1, 1)];
    end
else % half big half small
    Gn = [ratio * ones(Nh, 1); ones(Nh, 1)];
end

D = sqrt((L0^2 * sin(theta) * 0.0007) / (pi * N)); %L0^2 * sin(theta) = area of parallelogram 
D_all = Gn * D ; %particle size scaled so packing fraction is ~10^-3

xy_c_unit = [0; 1; 1; 0; 0; 0; 1; 1];
J = [1, cos(theta);
         0, sin(theta)];
    
xy_c_ori = reshape(xy_c_unit, Nvert, 2)'; 
xy_c_strain = J * xy_c_ori;

rng(pid); % random number generator with pid as the seed number

while true % generate random configuration with initial packing fraction 0.001, no particle is in contact with any other particle
    xy = rand(2 * N, 1) * L0;
    
    xy_c = reshape(xy_c_strain', Nvert * 2, 1) * L0;

    ol = overlap(xy, N, D_all, xy_c);
    if ol == 0
        break;
    end
end

xy_all = [xy; xy_c];
    
rsc = 1.001;
Pt_l = 1E-7; % low bound of total energy for the system to be deemed as jammed
Pt_h = 1.01 * Pt_l; % high bound of total energy for the system to be deemed as jammed

%Kp =  Kp_ratio * Kw;
Kpw = Kpw_ratio * Kp;
Kw = Kw_ratio * Kp;
dt_fire = 0.01 / sqrt(max([1, 4 * Kpw, 1/D_all(1)^2]));
Fthresh = 1E-14 * max([1, 4 * Kpw, 1/D_all(1)^2]);
Nt_fire = 1E7;
%%
success = 0;
D_l = -1.0;
D_h = -1.0;

% first compress the particles to a state close to jamming
while true
    if any(isnan(xy_all))
        break;
    end
    [xy_all, ~] = FIRE_Disk(xy_all, N, D_all, L0, Kp, Kpw, Fthresh, dt_fire, Nt_fire); % energy minimization, equivalently force balance
    [~, P] = StressTensor(xy_all, N, D_all, Kp, Kpw);
    if test == 1
        fprintf("P: %.5e   D: %.5e\n", P, D);
    end

    if P < Pt_l % unjammed state
        xy_old = xy_all;
        D_l = D;
        D = D * rsc;
        D_all = Gn * D;

        xy_all(1:N*2) = xy_all(1:N*2) * rsc; %I'm not sure if this part makes sense

    elseif P > Pt_h % over-compressed state
        D_h = D;
        D = (D_h + D_l) / 2.0;
        D_all = Gn * D;
        xy_all(1:N*2) = xy_old(1:N*2) * (D / D_h);
        
        if test == 1
            fprintf('D/D_l:  %.5e\n', D / D_h);
        end

        break;
    else % energy right within the threshold bounds
        success = 1;
        break;
    end
end
%% decrease particle growth(shrink) rate to faster find the jammed state
if success < 0.5
    rsc = 1.000001; 
    while true
        [xy_all, ~] = FIRE_Disk(xy_all, N, D_all, L0, Kp, Kpw, Fthresh, dt_fire, Nt_fire);
        [~, P] = StressTensor(xy_all, N, D_all, Kp, Kpw);
        if test == 1
            fprintf("P: %.5e   D: %.5e\n", P, D);
        end

        if P < Pt_l
            %fprintf('low\n');
            xy_old = xy_all;
            D_l = D;
            if D_h > 0.0
                %fprintf('D_l > 0\n')
                D = (D_h + D_l) / 2.0;
                D_h = -1.0;
            else
                D = D * rsc;
            end
            
            D_all = Gn * D;
            xy_all(1:N*2) = xy_old(1:N*2) * (D / D_l);

        elseif P > Pt_h
            D_h = D;
            D = (D_h + D_l) / 2.0;
            D_all = Gn * D;
            xy_all(1:N*2) = xy_old(1:N*2) * (D / D_h);
        
        else
            break;
        end
        if ((D_l > 0.0) && (abs(D_h / D_l - 1.0) < 1E-14))
            break;
        end
    end
end