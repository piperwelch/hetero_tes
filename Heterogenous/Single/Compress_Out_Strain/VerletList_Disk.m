function [VL, xy_save] = VerletList_Disk(xy, xy_save_old, N, D_all, first_call, VL_old)
%%
r_factor = 1.2;
r_cut = D_all(1);
r_list = r_factor * r_cut;
r_list_sq = r_list^2;
r_skin_sq = ((r_factor - 1) * r_cut)^2;
xy = xy(1:2*N);
%%
if first_call < 0.5
    dxy = xy - xy_save_old;
%     dxy = dxy - round(dxy / L) * L;
    dr = dxy(1:N).^2 + dxy(N+1:2*N).^2;
    if 4 * max(dr) < r_skin_sq
        VL = VL_old;
        xy_save = xy_save_old;
        return;
    end
end

%%
x = xy(1:N);
y = xy(N+1:2*N);
VL = zeros(N * 10, 2);
VL_counter = 0;

for n = 1:N-1
    for m = n+1:N
        dx = x(m) - x(n);
        if abs(dx) < r_list
            dy = y(m) - y(n);
            if abs(dy) < r_list                
                dr = dx^2 + dy^2;
                if dr < r_list_sq
                    VL_counter = VL_counter + 1;
                    VL(VL_counter, :) = [n, m];
                end
            end
        end
    end
end

VL = VL(1:VL_counter, :);

xy_save = xy;