function [Wlist, linklist, linklist_inner, ext_list] = ConstructTruss(Nvoxel_row)
%%
Nvoxel = Nvoxel_row * Nvoxel_row;
Nvert = (Nvoxel_row + 1)^2;

Wlist = zeros(4, Nvoxel);
count = 0;
for i = 1:Nvoxel_row
    for j = 1:Nvoxel_row
        count = count + 1;
        Wlist(:, count) = [i + j - 1; i + j; i + j + Nvoxel_row + 1; i + j + Nvoxel_row] + (i - 1) * Nvoxel_row;
    end
end

linklist = zeros(4 * Nvoxel, 2);
count = 0;
for i = 1:Nvoxel_row
    for j = 1:Nvoxel_row
        count = count + 2;
        idx1 = (i - 1) * (Nvoxel_row + 1) + j;
        linklist(count - 1, :) = [idx1, idx1 + 1];
        linklist(count, :) = [idx1, idx1 + Nvoxel_row + 1];
    end
    count = count + 1;
    linklist(count, :) = [idx1 + 1, idx1 + Nvoxel_row + 2];
end
for j = 1:Nvoxel_row
    count = count + 1;
    idx1 = Nvoxel_row * (Nvoxel_row + 1) + j;
    linklist(count, :) = [idx1, idx1 + 1];
end
linklist = linklist(1:count, :);

idx_outer = zeros(2 * Nvert, 1);
count = 0;
for j = 1:Nvoxel_row
    count = count + 2;
    idx = find(ismember(linklist, [j, j + 1], 'rows'));
    idx_outer(count - 1) = idx;
    idx = find(ismember(linklist, [j, j + 1] + Nvoxel_row * (Nvoxel_row + 1), 'rows'));
    idx_outer(count) = idx;
end
for i = 1:Nvoxel_row
    count = count + 2;
    idx1 = (i - 1) * (Nvoxel_row + 1) + 1;
    idx = find(ismember(linklist, [idx1, idx1 + Nvoxel_row + 1], 'rows'));
    idx_outer(count - 1) = idx;
    idx = find(ismember(linklist, [idx1, idx1 + Nvoxel_row + 1] + Nvoxel_row, 'rows'));
    idx_outer(count) = idx;
end
idx_outer = idx_outer(1:count);
linklist_inner = linklist;
linklist_inner(idx_outer, :) = [];

ext_list = [(1:Nvoxel_row+1)'; (Nvoxel_row+2:Nvoxel_row+1:(Nvoxel_row+1)^2)';...
            (Nvoxel_row+2:Nvoxel_row+1:(Nvoxel_row+1)^2)' + Nvoxel_row;...
            (Nvoxel_row*(Nvoxel_row+1)+2:(Nvoxel_row+1)^2-1)'];
ext_list = sort(ext_list);