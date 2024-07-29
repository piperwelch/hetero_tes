function [xy_p_all_comp, xy_c_all_comp] = compress(het_dir, xy_p_all, xy_c_all, N, Nvoxel_row, D_all, L0, L0_voxel, Kp, Kpw_ratio, Kw_ratio, test)

    cd(fullfile(het_dir, 'Compress'));
    n_pressure = 40; % number of pressures between 10^-7 and 10^-2 to sample
    tic
    [xy_p_all_comp, xy_c_all_comp, L0_voxel_comp, P_comp] = Compress(xy_p_all, xy_c_all, N, Nvoxel_row, D_all, L0, L0_voxel, Kp, Kpw_ratio, Kw_ratio, test);
    P_t_target = logspace(-7, -2, n_pressure)'; % target pressure values
    idx_target = zeros(n_pressure, 1);  % best pressure index
    for i = 1:n_pressure
        [~, idx] = min(abs(P_t_target(i) - P_comp));
        idx_target(i) = idx;
    end
    idx_target = unique(idx_target);
    P_comp = P_comp(idx_target);
    L0_voxel_comp = L0_voxel_comp(idx_target);
    xy_p_all_comp = xy_p_all_comp(:, idx_target');
    xy_c_all_comp = xy_c_all_comp(:, idx_target');
    toc
end