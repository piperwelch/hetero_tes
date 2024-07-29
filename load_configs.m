function [pos_unique, Kp, Kw_ratio, Kpw_ratio, Kw_ratio_tess, N, Gn] = load_configs(het_dir)
    dir_str = strcat(het_dir, '/Data/theta0.7854');
    pos_unique = readmatrix([dir_str,'/pos_N6.txt']);
    Gn = readmatrix([dir_str,'/Gn_N6.txt']);
    N = 6;

    % initialize variables & arrays
    Kp = 1;
    Kw_ratio = 0.1;
    Kpw_ratio = 1;
    Kw_ratio_tess = Kw_ratio * 0.5; 
end