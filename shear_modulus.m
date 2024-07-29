function [G_step, tess_Gdc] = shear_modulus(het_dir, xy_p_all_comp,xy_c_all_comp,  N, Nvoxel_row, D_all, L0, Kp, Kpw_ratio, Kw_ratio)

    tic
    cd(fullfile(het_dir, 'ShearModulus'));

    % ~20 sec to run 50 values of theta
    step = 1; % pressure index to measure at
    dgamma = 5e-9;

    Ntheta = 5;
    theta_all = linspace(0, pi / 2, Ntheta)';
    G_step = zeros(Ntheta, 1);

    for i = 1:Ntheta
        [sigma_shear, gamma_shear, G_all, xy_p_shear, xy_c_shear] = ...
                ShearModulus(xy_p_all_comp(:, step), xy_c_all_comp(:, step), N, Nvoxel_row, D_all, L0, Kp, Kpw_ratio, Kw_ratio, dgamma, theta_all(i));
        G_step(i) = G_all(1);
    end
    tess_info = FitSines(G_step', theta_all', 1); % note: inputs need to be N_trajectories x N_angles
    tess_Gamp = tess_info(4); % G_amp; as in G = G_amp * sin(4*th + th_0) + G_dc
    tess_Gdc = tess_info(5); % G_dc
    tess_th0 = tess_info(6); % theta_0 (in radians)
    toc
end