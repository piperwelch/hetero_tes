function isRattler = isRattler(x, y, xc, yc, D, N)
%%
ift = [2; 3; 4; 1];
lx = xc(ift) - xc;
ly = yc(ift) - yc;
L_all = sqrt(lx.^2 + ly.^2);
Cl_all = lx .* yc(ift) - ly .* xc(ift); % constant for line equation for four walls
%%
Z = zeros(N, 1);
for n = 1:N
    for w = 1:4
        % signed distance to Wall w, d > 0 if inside of the box
        if (lx(w) * y(n) - ly(w) * x(n) - Cl_all(w)) / L_all(w) < D(n) / 2
            Z(n) = Z(n) + 1;
        end
    end
    for m = n+1:N
        dx = x(m) - x(n);
        dy = y(m) - y(n);
        dr = sqrt(dx^2 + dy^2);
        if dr < 0.5 * (D(n) + D(m))
            Z(n) = Z(n) + 1;
            Z(m) = Z(m) + 1;
        end
    end
end

isRattler = (Z < 3);
%%
while true
    over = true;
    for n = 1:N
        if isRattler(n)
            continue;
        end
        
        theta = [];
        
        for w = 1:4
        % signed distance to Wall w, d > 0 if inside of the box
            if (lx(w) * y(n) - ly(w) * x(n) - Cl_all(w)) / L_all(w) < D(n) / 2
                theta = cat(1, theta, atan2(-lx(w) / L_all(w), ly(w) / L_all(w)));
            end
        end
        
        for m = 1:N
            if isRattler(m) || n == m
                continue;
            end
            dx = x(m) - x(n);
            dy = y(m) - y(n);
            dr = sqrt(dx^2 + dy^2);
            if dr < 0.5 * (D(n) + D(m))
                theta = cat(1, theta, atan2(dy, dx));
            end
        end
        
        Zc = size(theta, 1);
        if Zc < 3
            isRattler(n) = true;
            over = false;
            continue;
        else
            theta = sort(theta);
            if theta(Zc) - theta(1) < pi
                isRattler(n) = true;
                over = false;
            end
            if ~isRattler(n)
                for i = 1:Zc-1
                    if theta(i + 1) - theta(i) > pi
                        isRattler(n) = true;
                        over = false;
                        break;
                    end
                end
            end
        end
    end
    if over
        break;
    end
end 