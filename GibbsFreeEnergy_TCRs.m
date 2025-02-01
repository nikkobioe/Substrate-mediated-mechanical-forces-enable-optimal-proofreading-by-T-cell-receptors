%Gibbs Free energy calculations for two-pathway model and bell model for
%%TCR-pMHC binding

% General parameters
T_0 = 0.002; % TCR membrane density (nm^-2)
k_onT = 0.002./T_0; % normalized TCR on rate
F = linspace(0,100,1000); % transition state displacement force 0 pN - 100 pN
kT = 4.114; % Thermal energy = 4.114 pN*nm

%OT1 TCR - R4 MHC-1 slip bond
F_0sTs_OT1 = 10; % TCR-n-pMHC slip force (pN)
k_0sTs_OT1 = 2.0; % TCR-n-pMHC slip rate (s^-1)

k_offTs_OT1 = k_0sTs_OT1 .* exp(F ./ F_0sTs_OT1); % TCR-n-pMHC slip rate
K_D_OT1_R4 = k_offTs_OT1 ./ k_onT;
dG_OT1_R4 = kT .* log(K_D_OT1_R4);

%OT1 TCR - OVA MHC-I catch bond
k_0cTc_OT1 = 3.99; % TCR-c-pMHC catch rate (s^-1)
k_0sTc_OT1 = 0.41; % TCR-c-pMHC slip rate (s^-1)
F_0sTc_OT1 = 9.57; % TCR-c-pMHC slip force (pN)

k_offTc_OT1 = k_0cTc_OT1 .* exp(-F ./ F_0cTc_OT1) + k_0sTc_OT1 .* exp(F ./ F_0sTc_OT1); % catch bond off rate
K_D_OT1 = k_offTc_OT1 ./ k_onT; % dissociation constant
dG_OT1 = kT .* log(K_D_OT1); % Gibbs free energy

%OT1 TCR - A2 MHC-I catch bond
k_0cTc_A2 = 3.36; % TCR-c-pMHC catch rate (s^-1)
k_0sTc_A2 = 0.86; % TCR-c-pMHC slip rate (s^-1)
F_0cTc_A2 = 4.07; % TCR-c-pMHC catch force (pN)
F_0sTc_A2 = 11.8; % TCR-c-pMHC slip force (pN)

k_offTc_A2 = k_0cTc_A2 .* exp(-F ./ F_0cTc_A2) + k_0sTc_A2 .* exp(F ./ F_0sTc_A2); % catch bond off rate
K_D_A2 = k_offTc_A2 ./ k_onT; % dissociation constant
dG_A2 = kT .* log(K_D_A2); % Gibbs free energy

%OT1 TCR - G4 MHC-I catch bond
k_0cTc_G4 = 1.99; % TCR-c-pMHC catch rate (s^-1)
k_0sTc_G4 = 1.47; % TCR-c-pMHC slip rate (s^-1)
F_0cTc_G4 = 1.31; % TCR-c-pMHC catch force (pN)
F_0sTc_G4 = 3.96; % TCR-c-pMHC slip force (pN)

k_offTc_G4 = k_0cTc_G4 .* exp(-F ./ F_0cTc_G4) + k_0sTc_G4 .* exp(F ./ F_0sTc_G4); % catch bond off rate
K_D_G4 = k_offTc_G4 ./ k_onT; % dissociation constant
dG_G4 = kT .* log(K_D_G4); % Gibbs free energy

%OT1 TCR - E1 MHC-I catch bond
k_0cTc_E1 = 0.39; % TCR-c-pMHC catch rate (s^-1)
k_0sTc_E1 = 2.01; % TCR-c-pMHC slip rate (s^-1)
F_0cTc_E1 = 1.77; % TCR-c-pMHC catch force (pN)
F_0sTc_E1 = 4.33; % TCR-c-pMHC slip force (pN)

k_offTc_E1 = k_0cTc_E1 .* exp(-F ./ F_0cTc_E1) + k_0sTc_E1 .* exp(F ./ F_0sTc_E1); % catch bond off rate
K_D_E1 = k_offTc_E1 ./ k_onT; % dissociation constant
dG_E1 = kT .* log(K_D_E1); % Gibbs free energy

figure(1) % Gibbs free energy vs. Force
semilogx(F, dG_OT1./kT, F, dG_OT1_R4./kT, F, dG_A2./kT, F, dG_G4./kT, F, dG_E1./kT)
xlabel('Force (pN)')
ylabel('\Delta G (k_B T)')
xlim([min(F) max(F)]);
ylim([0 10]);
legend('OT1-TCR:R4-MHC-1', 'OT1-TCR:OVA-MHC-1', 'OT1-TCR:A2-MHC-1', 'OT1-TCR:G4-MHC-1', 'OT1-TCR:E1-MHC-1')

figure(2) % Bond lifetimes vs. Force
plot(F, 1./k_offTs_OT1, F, 1./k_offTc_OT1, F, 1./k_offTc_A2, F, 1./k_offTc_G4, F, 1./k_offTc_E1);
xlabel('Force (pN)')
ylabel('Bond lifetime (s)')
xlim([0 50]);
legend('OT1-TCR:R4-MHC-1', 'OT1-TCR:OVA-MHC-1', 'OT1-TCR:A2-MHC-1', 'OT1-TCR:G4-MHC-1', 'OT1-TCR:E1-MHC-1')

%Bond lifetimes

Tau_OVA = 1./k_offTc_OT1; 
Tau_R4 = 1./k_offTs_OT1;
Tau_A2 = 1./k_offTc_A2;
Tau_G4 = 1./k_offTc_G4;
Tau_E1 = 1./k_offTc_E1;

Tau_Matrix = [Tau_OVA; Tau_R4; Tau_A2; Tau_G4; Tau_G4; Tau_E1];
dG_Matrix_norm_kT = [dG_OT1./kT; dG_OT1_R4./kT; dG_A2./kT; dG_G4./kT; dG_E1./kT;];
