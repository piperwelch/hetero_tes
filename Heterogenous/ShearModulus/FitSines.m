function [sin_info] = FitSines(data_arr, theta_list, visu)

%{ 
----------- inputs ------------
data_arr: (n_data x n_angles), each row is a unique g(theta) or lambda(theta), where there are n_angles value of theta
theta_list: (1 x n_angles): list of angles theta
visu: (0 or 1 bool) visualize fits (NOTE: need to press key to proceed in for loop)

----------- outputs -----------
% sin_info: (n_data x 6): [init a, init dc, init th0, fitted a, fitted dc, fitted th0]
defined as G(theta) = G_a * sin(4*theta + theta_0) + G_dc
%}

sin_info = zeros(size(data_arr, 1), 6); % store values of {amp/dc/theta0} for each row

for kk = [1:size(data_arr, 1)] % loop through elements in the array
    current_traj = squeeze(data_arr(kk, :)); % current f(theta) trace
    current_amp = (max(current_traj, [], 2) - min(current_traj, [], 2))/2; % estimated amplitude
    current_dc = ((max(current_traj, [], 2) - current_amp) + (min(current_traj, [], 2) + current_amp)) / 2; % estimated dc

    [~, max_ind] = min(abs(current_traj));
    current_thzero = theta_list(max_ind); % estimated theta

    if current_thzero <= 0.373737*pi % weird ansatz wrapping fix for theta_0
        current_thzero = (current_thzero-3*pi/8)*4;
    else
        current_thzero = (current_thzero-7*pi/8)*4;
    end

    % store guesses for amp/dc/theta
    sin_info(kk, 1:3) = [current_amp current_dc current_thzero];
    
    % fit constants using matlab's "fit" function, if error, continue and store NaN
    ft = fittype('sinfit(x, a, dc, th0)');
    try
        f = fit(theta_list', current_traj', ft, 'startpoint', [current_amp, current_dc, current_thzero], 'lower', [0, -Inf, 0], 'upper', [Inf, Inf, 2*pi]);
    catch
        sin_info(kk, 4:6) = [NaN NaN NaN];
        fprintf("fitting error on %.0f \n", kk)            
        continue
    end

    fit_th0 = f.th0;        
    sin_info(kk, 4:6) = [f.a f.dc fit_th0];

    %if visu % plot visualization
    %    figure(3), clf, hold on        
    %    str_type = 'ko-';
    %
    %    plot(theta_list ./ pi, current_traj, str_type, 'linewidth', 0.3)
    %    plot(theta_list ./ pi, sinfit(theta_list, f.a, f.dc, fit_th0), 'r--', 'linewidth', 1)
    %    plot(theta_list ./ pi, sinfit(theta_list, current_amp, current_dc, current_thzero), 'c--', 'linewidth', 1)
    %    legend("data", "fit", "init guess")
    %    xlabel('\theta / \pi')
    %    ylabel('\lambda_i')
    %    set(gca, 'fontsize', 14)
    %    title(num2str(kk))
%             pause(0.5)
    %end



end
end
