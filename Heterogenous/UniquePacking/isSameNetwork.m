function is_same = isSameNetwork(pos_all_1, pos_all_2, N)
%%
is_same = false;
Gn = [1.4 * ones(N / 2, 1); ones(N / 2, 1)];
%%
D1 = pos_all_1(1) * Gn;
x1 = pos_all_1(2:N+1);
y1 = pos_all_1(N+2:2*N+1);
xc1 = pos_all_1(2*N+2:2*N+5);
yc1 = pos_all_1(2*N+6:2*N+9);
D2 = pos_all_2(1) * Gn;
x2 = pos_all_2(2:N+1);
y2 = pos_all_2(N+2:2*N+1);
xc2 = pos_all_2(2*N+2:2*N+5);
yc2 = pos_all_2(2*N+6:2*N+9);
[Nnr1, Nnr_b_1, Nnr_s_1, Z_all1, A1] = ContactNetwork_forCheck(x1, y1, xc1, yc1, D1, N);
[Nnr2, Nnr_b_2, Nnr_s_2, Z_all2, A2] = ContactNetwork_forCheck(x2, y2, xc2, yc2, D2, N);
%%
if Nnr1 ~= Nnr2 || Nnr_b_1 ~= Nnr_b_2 || Nnr_s_1 ~= Nnr_s_2 || Z_all1 ~= Z_all2
    return;
end

if abs(det(A1) - det(A2)) > 1e-10
    return;
end

eigD1 = eig(A1);
eigD2 = eig(A2);
if sum((eigD1 - eigD2).^2) > 1e-10
    return;
end

if N > 10
    is_same = true;
    return;
end
%%
A1 = A1(1:Nnr1, 1:Nnr1);
A2 = A2(1:Nnr2, 1:Nnr2);

Nnr = Nnr1;
Nnr_b = Nnr_b_1;
Nnr_s = Nnr_s_1;
I = eye(Nnr);
pr_b = perms(1:Nnr_b);
pr_s = perms(1:Nnr_s) + Nnr_b;
perm_num = size(pr_b, 1) * size(pr_s, 1);
pr = zeros(perm_num, Nnr);
pr_count = 0;
for i = 1:size(pr_b, 1)
    for j = 1:size(pr_s, 1)
        pr_count = pr_count + 1;
        pr(pr_count, :) = [pr_b(i, :), pr_s(j, :)];
    end
end

for i = 1:perm_num
    P = I(pr(i, :), :);
    if isequal(P * A1, A2 * P)
        is_same = true;
        return;
    end
end

% if gamma > 0
%     for i = 1:perm_num
%         pp = I(pr(i, :), :);
%         P = eye(Nnr + 4);
%         P(1:Nnr, 1:Nnr) = pp;
%         if isequal(P * A1, A2 * P)
%             is_same = true;
%             return;
%         end
%     end
% else
%     for i = 1:perm_num
%         pp = I(pr(i, :), :);
%         P1 = eye(Nnr + 4);
%         P1(1:Nnr, 1:Nnr) = pp;
%         P2 = P1;
%         P2(Nnr + 2, Nnr + 2) = 0;
%         P2(Nnr + 2, Nnr + 4) = 1;
%         P2(Nnr + 4, Nnr + 2) = 1;
%         P2(Nnr + 4, Nnr + 4) = 0;
%         if isequal(P1 * A1, A2 * P1) || isequal(P2 * A1, A2 * P2)
%             is_same = true;
%             return;
%         end
%     end
% end
