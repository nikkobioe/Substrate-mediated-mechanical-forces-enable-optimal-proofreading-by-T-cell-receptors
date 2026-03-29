%% Substrate mediated mechanical forces enable optimal kinetic proofreading by T-cell receptors 

%% Jeffreys & Shankar et al.
%% TCR and LFA-1 forces steady-state analysis - solve roots

% Parameters
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
gamma = 3e-2; % T-cell membrane surface tension (pN/nm)
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

% LFA-1/TCR forces and APC surface displacement as a function of APC stiffness
x0 = [12.5, 10.5]; % set initial guesses for F_ic and F_Tc

options = optimoptions('lsqnonlin', 'MaxFunctionEvaluations', 5000, 'MaxIterations', 5000);
sigma_apc = logspace(-3, 3, 1000); % Substrate stiffness
F_ic = zeros(1,1000); % Initialize array to store results for LFA-1 forces
F_Tc = zeros(1,1000); % Initialize array to store results for TCR-c-pMHC forces
F_mc = zeros(1,1000); % Initialize array to store results for T-cell membrane forces
u = zeros(1,1000); % Initialize array to store results for substrate displacement
dh = zeros(1,1000); % Initialize array to store results for membrane displacement
ic = zeros(1,1000); % Initialize array to store results for LFA-1-ICAM-1 bond density
v_a = zeros(1,1000); % Initialize array to store results for actin retrograde flow velocity
PI_c = zeros(1,1000); % Initialize array to store results for TCR-c-pMHC bond density
yc = zeros(1,1000); % Initalize array to store results for average LFA-1-ICAM-1 bond stretch 
k_offic = zeros(1,1000); % Initialize array to store results for LFA-1 off rate
k_offTc = zeros(1,1000); % Initialize array to store results for TCR off rate

for i = 1:1000
    % Define the objective function for lsqnonlin
    objective_fun = @(x) TCRcatchv2(sigma_apc(i), x(1), x(2), k_0ci, F_0ci, k_0si, F_0si, k_0cTc, F_0cTc, k_0sTc, F_0sTc, k_onT, k_oni, i_0, v_0, ...
    sigma_a, p_c, T_0, gamma, k_T, k_i, l_0i, l_0T, R, a_c, h_0);
    
    % Solve using lsqnonlin, starting with x0
    x = lsqnonlin(objective_fun, x0, [], [], options); 
    F_ic(i) = x(1);
    F_Tc(i) = x(2);
    
    % Calculate adhesion kinetic, proofreading, and substrate mechanics based on the obtained F_ic and F_Tc
    % Calculate LFA-1 adhesion kinetics and dissociation constants
    k_offic(i) = k_0ci .* exp(-F_ic(i) ./ F_0ci) + k_0si .* exp(F_ic(i) ./ F_0si); % LFA-1 catch rate
    k_offTc(i) = k_0cTc .* exp(-F_Tc(i) ./ F_0cTc) + k_0sTc .* exp(F_Tc(i) ./ F_0sTc); % TCR-c-pMHC catch rate
    K_dc = k_offTc(i) ./ k_onT; % TCR-c-pMHC dissociation constant
    K_dic = k_offic(i) ./ k_oni; % LFA-1 dissociation constant cognate
    ic(i) = i_0 ./ (1 + K_dic); % LFA-1-ICAM-1 bond density
    v_a(i) = v_0 .* (1 - ic(i) .* F_ic(i) ./ sigma_a); % actin retrograde flow velocity
    yc(i) = v_a(i) ./ k_offic(i); % mean LFA-1-ICAM-1 bond stretch length

    % TCR-ligand bond density
    PI_c(i) = 0.5 .* (p_c + T_0 + K_dc) - 0.5 .* sqrt((p_c + T_0 + K_dc).^2 - 4 * p_c * T_0);

    % Solve for T-cell membrane displacement
    k_m = gamma / R^2; % T-cell membrane force density (pN/nm^3)
    F_mc(i) = PI_c(i) .* F_Tc(i); % membrane generated stresses (pN/nm^2)
    dh(i) = F_mc(i) ./ k_m; % membrane displacement (nm)
    
    % Solve for APC membrane displacement
    u(i) = (a_c ./ sigma_apc(i)) .* (PI_c(i) .* F_Tc(i) + ic(i) .* F_ic(i));
end

% Plot results
figure(1);
semilogx(sigma_apc, F_ic, '-o');
title('LFA-1 force');
ylabel('LFA-1 force (pN)');
xlabel('APC stiffness (pN/nm)');
xlim([1e-3 1e3]);

figure(2);
semilogx(sigma_apc, u, '-o');
title('APC substrate displacement');
ylabel('Displacement (nm)');
xlabel('APC stiffness (pN/nm)');
xlim([1e-3 1e3]);

figure(3);
semilogx(sigma_apc, ic .* 1e6, '-o');
title('LFA-1-ICAM-1 bond density');
ylabel('Bond density (um^-^2)');
xlabel('APC stiffness (pN/nm)');
xlim([1e-3 1e3]);

figure(4);
semilogx(sigma_apc, 1 ./ k_offic, '-o');
title('LFA-1-ICAM-1 bond lifetime');
ylabel('Bond lifetime (s)');
xlabel('APC stiffness (pN/nm)');
xlim([1e-3 1e3]);

figure(5);
semilogx(sigma_apc, F_Tc, '-o');
title('TCR-c-pMHC force');
ylabel('TCR force (pN)');
xlabel('APC stiffness (pN/nm)');
xlim([1e-3 1e3]);

figure(6);
semilogx(sigma_apc, PI_c .* 1e6, '-o');
title('TCR-c-pMHC bond density');
ylabel('Bond density (um^-^2)');
xlabel('APC stiffness (pN/nm)');
xlim([1e-3 1e3]);

figure(7);
semilogx(sigma_apc, 1 ./ k_offTc, '-o');
title('TCR-c-pMHC bond lifetime');
ylabel('Bond lifetime (s)');
xlabel('APC stiffness (pN/nm)');
xlim([1e-3 1e3]);

figure(8);
semilogx(sigma_apc, F_mc, '-o');
title('T-cell membrane stress');
ylabel('Membrane stress (pN/nm^2)');
xlabel('APC stiffness (pN/nm)');
xlim([1e-3 1e3]);

figure(9);
semilogx(sigma_apc, dh, '-o');
title('T-cell membrane displacement');
ylabel('Displacement (nm)');
xlabel('APC stiffness (pN/nm)');
xlim([1e-3 1e3]);

figure(10);
semilogx(sigma_apc, v_a, '-o');
title('actin retrograde flow');
ylabel('velocity (nm/s)');
xlabel('APC stiffness (pN/nm)');
xlim([1e-3 1e3]);

% kinetic proofreading with limited signaling
alpha_c = k_p ./ (k_p + k_offTc); % phosphorylation rate competition
beta_c = k_offTc ./ (k_offTc + phi); % signal decay rate competition
R_c = zeros(theta, 1000); % Initialize array for signaling catch bond TCRs

figure(11)
hold on;
for i = 1:theta
    R_c(i,:) = alpha_c.^(i) .* beta_c .* PI_c;
    semilogx(sigma_apc, R_c(i,:).*1e6, '-o');
end
title('Signaling TCR-c-pMHC bond density');
ylabel('Signaling complex density (um^-^2)');
xlabel('APC stiffness (pN/nm)');
xlim([1e-3 1e3]);
set(gca, 'XScale', 'log'); % enforce logscale
hold off;

% additional arrays for results
tau_Tc = 1 ./ k_offTc; % TCR-c-pMHC bond lifetime
tau_ic = 1 ./ k_offic; % LFA-1-ICAM-1 bond lifetime
PI_c_permicronsquared = PI_c.*1e6;
R_c_permicronsquared = R_c.*1e6;
ic_permicronsquared = ic.*1e6;
K_Dc = k_offTc./k_onT;
K_Dic = k_offic./k_oni;
