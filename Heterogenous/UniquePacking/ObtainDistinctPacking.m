function [pos_unique_act, pos_unique, config_idx, appear] = ObtainDistinctPacking(N, pos_all)

%{
Generates jammed configs in boundaries at fixed theta (rhombus with length L)

input:
------
N (sc): particle number
pos_all (Nconfigs x 2N+9): first input is diameter, then particle postions, then box positions
-----

output:
------
pos_unique_act (N_unique x 2N+9): unique rows of pos_all according to contact network and phi (organized also as diameter then particle positions then box positions)
pos_unique: intermediate result, only unique configs. based on phi (packing volume)
config_idx (Nconfigs x 1): which unique config class each initial input belongs to
appear: loosely, if contact networks are same
------
%}

%%
D_all = pos_all(:, 1);
[D_unique, ~, d_ia] = unique(D_all);
D_num = size(D_unique, 1);
fprintf('d_ia:  %d\n', d_ia);
fprintf('D_num: %d\n', D_num);
pos_unique = [];
unique_count = 0;
config_idx = zeros(size(pos_all, 1), 1);

%%
for it = 1:D_num
    fprintf('it: %d  %.4f\n', it, it / D_num);
    % idx = find(abs(L_all - L_unique(it)) < 1e-8 / N);
    idx = find(d_ia == it);
    fprintf('idx(1):  %d\n', idx(1));
    unique_count = unique_count + 1;
    config_idx(idx(1)) = unique_count;
    a = pos_all(idx(1), :)';
    pos_unique = cat(1, pos_unique, a');
    pnum = length(idx) - 1;
    if pnum == 0
        continue;
    end
    idx(1) = [];
    while true % check if a configuration with the same L is a different configuration
        is_same_phi = (zeros(pnum, 1) > 1);
        for i = 1:pnum
            is_same_phi(i) = isSameNetwork(a, pos_all(idx(i), :)', N);
        end
        config_idx(idx(is_same_phi)) = unique_count;
        idx = idx(~is_same_phi);
        if size(idx, 1) == 0
            break;
        else
            pnum = length(idx) - 1;
            unique_count = unique_count + 1;
            config_idx(idx(1)) = unique_count;
            a = pos_all(idx(1), :)';
            pos_unique = cat(1, pos_unique, a');
            if pnum == 0
                break;
            end
            idx(1) = [];
        end
    end
end

%%
is_diff = zeros(unique_count * (unique_count - 1) / 2, 1);
count = 0;
unique_count_act = unique_count;
appear = zeros(unique_count, 1);
for i = 1:unique_count-1
    if appear(i) == 1
        continue;
    end
    a = pos_unique(i, :)';
    for j = i+1:unique_count
        count = count + 1;
        if appear(j) == 1
            continue;
        end
        is_same = isSameNetwork(a, pos_unique(j, :)', N);
        if ~is_same
            is_diff(count) = 1;
        else
            unique_count_act = unique_count_act - 1;
            appear(j) = 1;
            config_idx(config_idx == j) = i;
        end
    end
end
pos_unique_act = pos_unique((appear < 0.5), :);