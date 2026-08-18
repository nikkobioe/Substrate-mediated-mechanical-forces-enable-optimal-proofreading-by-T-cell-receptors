%% Parameter Sensitivity Analysis (Based on Bangasser et al. 2013)
clear; clc; close all;

% 1. Define Baseline Parameters in a Struct (29 Parameters)
param_names = {'T_0','i_0','p_c','p_n','k_T','k_i','k_oni','k_0ci','k_0si', ...
               'F_0ci','F_0si','k_onT','k_0cTc','k_0sTc','F_0cTc','F_0sTc', ...
               'k_0sTs','F_0sTs','gamma','sigma_a','v_0','h_0','l_0T', ...
               'l_0i','R','a_c','theta','k_p','phi'};
num_params = length(param_names);

params.T_0 = 0.002;
params.i_0 = 0.0004;
params.p_c = 0.002;
params.p_n = 0.02;
params.k_T = 0.6;
params.k_i = 0.2;
params.k_oni = 0.1;
params.k_0ci = 2.5;
params.k_0si = 1.4e-2;
params.F_0ci = 6.2;
params.F_0si = 3.9;
params.k_onT = 0.01 / params.T_0;
params.k_0cTc = 3.99;
params.k_0sTc = 0.41;
params.F_0cTc = 3.03;
params.F_0sTc = 9.57;
params.k_0sTs = 2.0;
params.F_0sTs = 10;
params.gamma = 3e-2;
params.sigma_a = 0.001;
params.v_0 = 60;
params.h_0 = 55;
params.l_0T = 15;
params.l_0i = 40;
params.R = 50;
params.a_c = 3e4;
params.theta = 10;
params.k_p = 1;
params.phi = 0.1;

% Define APC Stiffness Evaluation Array
sigma_apc = logspace(-3, 3, 100);

% Define the fold-change multipliers (c) as per the paper
c_vals = [0.25, 0.5, 1, 2, 4];
log_c = log10(c_vals);

% Initialize 29x29 Sensitivity Matrices (S values)
S_stiffness_matrix = zeros(num_params, num_params);
S_precision_matrix = zeros(num_params, num_params);

disp('Running Single and Pairwise Sensitivity Analysis...');

for i = 1:num_params
    for j = 1:num_params
        
        opt_stiffness_array = zeros(1, length(c_vals));
        max_precision_array = zeros(1, length(c_vals));
        
        for k = 1:length(c_vals)
            c = c_vals(k);
            temp_params = params;
            
            pname_i = param_names{i};
            pname_j = param_names{j};
            
            if i == j
                % SINGLE PARAMETER (Diagonal)
                temp_params.(pname_i) = params.(pname_i) * c;
            elseif j > i
                % DUAL PARAMETER - Same Direction (Upper Triangle)
                temp_params.(pname_i) = params.(pname_i) * c;
                temp_params.(pname_j) = params.(pname_j) * c;
            else
                % DUAL PARAMETER - Opposite Direction (Lower Triangle)
                temp_params.(pname_i) = params.(pname_i) * c;
                temp_params.(pname_j) = params.(pname_j) / c;
            end
            
            % Evaluate model
            precision_curve = evaluate_precision(temp_params, sigma_apc);
            
            % Find peak precision and associated optimal stiffness
            [max_val, max_idx] = max(precision_curve);
            
            max_precision_array(k) = max_val;
            opt_stiffness_array(k) = sigma_apc(max_idx);
        end
        
        % Calculate S (log-log slope) using polyfit
        % 1. Sensitivity of Optimal Stiffness
        p_stiff = polyfit(log_c, log10(opt_stiffness_array), 1);
        S_stiffness_matrix(i, j) = p_stiff(1);
        
        % 2. Sensitivity of Maximum Precision
        p_prec = polyfit(log_c, log10(max_precision_array + eps), 1);
        S_precision_matrix(i, j) = p_prec(1);
    end
    fprintf('Completed Parameter Row %d of 29 (%s)\n', i, param_names{i});
end

%% Plotting the Results
figure('Name', 'Sensitivity of Optimal Substrate Stiffness (S)', 'Position', [100, 100, 900, 800]);
h1 = heatmap(param_names, param_names, S_stiffness_matrix);
h1.Title = 'Sensitivity (S) of Optimal Stiffness to Dual Parameter Changes';
h1.XLabel = 'Parameter j';
h1.YLabel = 'Parameter i';
h1.Colormap = redbluecmap; % Custom or built-in diverging colormap for pos/neg slopes

figure('Name', 'Sensitivity of Maximum Precision (S)', 'Position', [150, 150, 900, 800]);
h2 = heatmap(param_names, param_names, S_precision_matrix);
h2.Title = 'Sensitivity (S) of Max Precision to Dual Parameter Changes';
h2.XLabel = 'Parameter j';
h2.YLabel = 'Parameter i';
h2.Colormap = redbluecmap; 

% Export 
writematrix(S_stiffness_matrix, 'S_optimal_stiffness_29x29.csv');
writematrix(S_precision_matrix, 'S_max_precision_29x29.csv');

%% Helper Function
function precision = evaluate_precision(p, sigma_apc)
    options = optimoptions('lsqnonlin', 'Display', 'off', ...
                           'MaxFunctionEvaluations', 500, 'MaxIterations', 500);
    num_points = length(sigma_apc);
    R_c = zeros(1, num_points);
    R_s = zeros(1, num_points);
    
    x0_catch = [12.5, 10.5];
    x0_slip = [12.5, 0];
    
    for k = 1:num_points
        % Catch bond
        obj_catch = @(x) TCRcatchv2(sigma_apc(k), x(1), x(2), p.k_0ci, p.F_0ci, p.k_0si, p.F_0si, ...
            p.k_0cTc, p.F_0cTc, p.k_0sTc, p.F_0sTc, p.k_onT, p.k_oni, p.i_0, p.v_0, ...
            p.sigma_a, p.p_c, p.T_0, p.gamma, p.k_T, p.k_i, p.l_0i, p.l_0T, p.R, p.a_c, p.h_0);
        x_catch = lsqnonlin(obj_catch, x0_catch, [], [], options);
        F_Tc = x_catch(2);
        
        k_offTc = p.k_0cTc * exp(-F_Tc / p.F_0cTc) + p.k_0sTc * exp(F_Tc / p.F_0sTc);
        K_dc = k_offTc / p.k_onT;
        PI_c = 0.5 * (p.p_c + p.T_0 + K_dc) - 0.5 * sqrt((p.p_c + p.T_0 + K_dc)^2 - 4 * p.p_c * p.T_0);
        alpha_c = p.k_p / (p.k_p + k_offTc);
        beta_c = k_offTc / (k_offTc + p.phi);
        R_c(k) = (alpha_c^p.theta) * beta_c * PI_c; 
        
        % Slip bond
        obj_slip = @(x) TCRslipv2(sigma_apc(k), x(1), x(2), p.k_0ci, p.F_0ci, p.k_0si, p.F_0si, ...
            p.k_0sTs, p.F_0sTs, p.k_onT, p.k_oni, p.i_0, p.v_0, ...
            p.sigma_a, p.p_n, p.T_0, p.gamma, p.k_T, p.k_i, p.l_0i, p.l_0T, p.R, p.a_c, p.h_0);
        x_slip = lsqnonlin(obj_slip, x0_slip, [], [], options);
        F_Ts = x_slip(2);
        
        k_offTs = p.k_0sTs * exp(F_Ts / p.F_0sTs);
        K_ds = k_offTs / p.k_onT;
        PI_s = 0.5 * (p.p_n + p.T_0 + K_ds) - 0.5 * sqrt((p.p_n + p.T_0 + K_ds)^2 - 4 * p.p_n * p.T_0);
        alpha_s = p.k_p / (p.k_p + k_offTs);
        beta_s = k_offTs / (k_offTs + p.phi);
        R_s(k) = (alpha_s^p.theta) * beta_s * PI_s;
    end
    precision = R_c ./ (R_s + eps);
end