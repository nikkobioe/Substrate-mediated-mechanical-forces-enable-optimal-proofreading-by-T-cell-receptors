%% Phase diagram of proofreading precision
% Evaluates precision at theta = 10 over a grid of APC stiffness and membrane tension

% Shared Parameters
T_0 = 0.002; % TCR membrane density (nm^-2)
i_0 = 0.0004; % LFA-1 membrane density (nm^-2)
p_c = 0.002; % c-pMHC surface density (nm^-2)
p_n = 0.02; % n-pHC surface density (nm^-2)
k_T = 0.6; % TCR-pMHC bond stiffness (pN/nm)
k_i = 0.2; % LFA-1-ICAM-1 bond stiffness (pN/nm)
k_oni = 0.1; % LFA-1 on rate (s^-1)
k_0ci = 2.5; % LFA-1 catch rate (s^-1)
k_0si = 1.4e-2; % LFA-1 slip rate (s^-1)
F_0ci = 6.2; % LFA-1 catch force (pN)
F_0si = 3.9; % LFA-1 slip force (pN)
k_onT = 0.01/T_0; % TCR on rate per bond(nm^2*s^-1)
k_0cTc = 3.99; % TCR-c-pMHC catch rate (s^-1)
k_0sTc = 0.41; % TCR-c-pMHC slip rate (s^-1)
F_0cTc = 3.03; % TCR-c-pMHC catch force (pN)
F_0sTc = 9.57; % TCR-c-pMHC slip force (pN)
k_0sTs = 2.0; % TCR-n-pMHC slip rate (s^-1)
F_0sTs = 10; % TCR-n-pMHC slip force (pN)
sigma_a = 0.001; % actomyosin active stress (pN/nm^2)
v_0 = 60; % actomyosin unloaded retrograde flow velocity (nm/s)
h_0 = 55; % equilibrium separation distance (CD45, nm)
l_0T = 15; % TCR-pMHC complex bond length (nm)
l_0i = 40; % LFA-1-ICAM-1 complex bond length (nm)
R = 50; % TCR microcluster radii (nm)
a_c = 3e4; % contact area (nm^2)
theta = 10; % no. of proofreading sequences
k_p = 1; % tyrosine phosphorylation rate (s^-1)
phi = 0.1; % TCR-pMHC signaling complex decay rate

% Grid Setup for Phase Diagram
N = 100; % Grid resolution 
sigma_apc_vals = logspace(log10(0.001), log10(1000), N); % APC stiffness (pN/nm)
gamma_vals = logspace(log10(0.001), log10(10), N); % T-cell membrane tension (pN/nm)

precision_matrix = zeros(N, N);

% Suppress lsqnonlin output for cleaner execution
options = optimoptions('lsqnonlin', 'Display', 'off', 'MaxFunctionEvaluations', 1000, 'MaxIterations', 1000);

% Initial guesses based on previous 1D parameter sweeps
x0_catch = [12.5, 10.5]; 
x0_slip = [12.5, 0]; 

disp('Calculating phase diagram. This may take a few moments...');

for i = 1:N
    for j = 1:N
        sigma_apc = sigma_apc_vals(i);
        gamma = gamma_vals(j);
        
        %% --- Cognate Catch Bond Model ---
        obj_catch = @(x) TCRcatchv2(sigma_apc, x(1), x(2), k_0ci, F_0ci, k_0si, F_0si, ...
            k_0cTc, F_0cTc, k_0sTc, F_0sTc, k_onT, k_oni, i_0, v_0, ...
            sigma_a, p_c, T_0, gamma, k_T, k_i, l_0i, l_0T, R, a_c, h_0);
        
        x_catch = lsqnonlin(obj_catch, x0_catch, [], [], options);
        F_Tc = x_catch(2); % Extract TCR-c-pMHC force
        
        % Calculate cognate proofreading metrics
        k_offTc = k_0cTc * exp(-F_Tc / F_0cTc) + k_0sTc * exp(F_Tc / F_0sTc);
        K_dc = k_offTc / k_onT;
        PI_c = 0.5 * (p_c + T_0 + K_dc) - 0.5 * sqrt((p_c + T_0 + K_dc)^2 - 4 * p_c * T_0);
        
        alpha_c = k_p / (k_p + k_offTc);
        beta_c = k_offTc / (k_offTc + phi);
        R_c = (alpha_c^theta) * beta_c * PI_c; % Signaling TCR-c-pMHC bond density
        
        %% --- Non-Cognate Slip Bond Model ---
        obj_slip = @(x) TCRslipv2(sigma_apc, x(1), x(2), k_0ci, F_0ci, k_0si, F_0si, ...
            k_0sTs, F_0sTs, k_onT, k_oni, i_0, v_0, ...
            sigma_a, p_n, T_0, gamma, k_T, k_i, l_0i, l_0T, R, a_c, h_0);
        
        x_slip = lsqnonlin(obj_slip, x0_slip, [], [], options);
        F_Ts = x_slip(2); % Extract TCR-n-pMHC force
        
        % Calculate non-cognate proofreading metrics
        k_offTs = k_0sTs * exp(F_Ts / F_0sTs);
        K_ds = k_offTs / k_onT;
        PI_s = 0.5 * (p_n + T_0 + K_ds) - 0.5 * sqrt((p_n + T_0 + K_ds)^2 - 4 * p_n * T_0);
        
        alpha_s = k_p / (k_p + k_offTs);
        beta_s = k_offTs / (k_offTs + phi);
        R_s = (alpha_s^theta) * beta_s * PI_s; % Signaling TCR-n-pMHC bond density
        
        %% --- Precision Calculation ---
        % Defined as ratio of cognate to non-cognate signaling density
        % max() prevents division by zero and eliminates NaN values
        precision_matrix(j, i) = R_c / max(R_s, 1e-12);
    end
end

%% Plotting the Phase Diagram
[X, Y] = meshgrid(sigma_apc_vals, gamma_vals);

figure('Name', 'Precision Phase Diagram');
contourf(X, Y, log10(precision_matrix), 30, 'LineStyle', 'none');
c = colorbar;
c.Label.String = 'Log_{10}(Precision)';
xlabel('APC Surface Stiffness, \sigma_{APC} (pN/nm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('T-cell Membrane Tension, \gamma (pN/nm)', 'FontSize', 12, 'FontWeight', 'bold');
title('Proofreading Precision (\theta = 10)', 'FontSize', 14);
colormap('turbo'); % 'parula' 'cividis' or 'viridis' are also excellent alternatives
set(gca, 'XScale', 'log', 'YScale', 'log'); % enforce logscale for both axes
xlim([0.001 1000]);
ylim([0.001 10]);