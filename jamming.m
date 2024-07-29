function [xy_p_all, xy_c_all, D_all, Wlist, linklist, linklist_inner, ext_list, L0_voxel] = ...
    jamming(het_dir, conid, pos_unique, N, Gn, Kp, Kw_ratio, Kpw_ratio, Kw_ratio_tess)
    
    cd(fullfile(het_dir, 'Jamming'));
    Nvoxel_row = 4;
    Nvoxel = Nvoxel_row * Nvoxel_row;
    xy_voxels = zeros(size(pos_unique, 1), 2*N+8);
    D_every = pos_unique(:, 1);
    
    %rotations single voxels to all have bottom edge as horizontal (mostly formality)
    for j = 1:size(pos_unique, 1)
            xy_voxels(j, :) = RotateConfig_Voxel(N, pos_unique(j, 2:2*N+9), 1, 4,  [1; 2; 3; 4]);
    end
    
    %xy_p_all = 1; xy_c_all = 1; D_all= 1; Wlist = 1; linklist = 1; linklist_inner = 1; ext_list = 1; 
    % % create conid array - elements are index of config you want in that spot in array
    tessname = '11_config1';
    %conid = ones(Nvoxel, 1);
    %conid(1) = 1; 
    %conid(2) = 2; conid(3) = 5; conid(4) = 4; conid(5) = 5; conid(6) = 6; 
    %conid(7) = 7; conid(8) = 8; conid(9) = 9; conid(10) = 10; conid(11) = 11; 
    %conid(12) = 12; conid(13) = 13; conid(14) = 14; 
    % % tessellate system -> input conid, jamming uses it to place single voxels into array
    test = 0; L0 = 1; 
    [Wlist, linklist, linklist_inner, ext_list] = ConstructTruss(Nvoxel_row);
    [xy_p_all, xy_c_all, L0_voxel, D_all] = Jamming(conid, xy_voxels, N, D_every, Gn, L0, Nvoxel_row, Kp, Kw_ratio_tess, Kpw_ratio, test);
    % 
    % % plot configuration to check: 
    % PlotConfig(xy_p_all, xy_c_all, N, Nvoxel_row, D_all, Wlist, linklist, ext_list, 1)

end