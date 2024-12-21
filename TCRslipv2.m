function F = TCRslipv2(sigma_apc, F_ic, F_Tn, k_0ci, F_0ci, k_0si, F_0si, k_0sTs, F_0sTs, k_onT, k_oni, i_0, v_0, ...
    sigma_a, p_n, T_0, gamma, k_T, k_i, l_0i, l_0T, R, a_c, h_0)
    
    % Calculate LFA-1 adhesion kinetics and dissociation constants
    k_offic = k_0ci * exp(-F_ic / F_0ci) + k_0si * exp(F_ic / F_0si);
    k_offTs = k_0sTs * exp(F_Tn / F_0sTs);
    K_ds = k_offTs / k_onT;
    K_dic = k_offic / k_oni;
    ic = i_0 / (1 + K_dic);
    v_a = v_0 * (1 - ic * F_ic / sigma_a);
    yc = v_a / k_offic;

    % TCR-ligand bond density
    PI_n = 0.5 * (p_n + T_0 + K_ds) - 0.5 * sqrt((p_n + T_0 + K_ds)^2 - 4 * p_n * T_0);

    % Solve for T-cell membrane displacement
    k_m = gamma / R^2; % T-cell membrane force density (pN/nm^3)
    F_mc = PI_n * F_Tn; % membrane generated stresses (pN/nm^2)
    dh = F_mc / k_m; % membrane displacement (nm)

    % Solve for APC membrane displacement
    u = (a_c / sigma_apc) * (PI_n * F_Tn + ic * F_ic);

    % Calculate LFA-1 and TCR forces
    F1 = k_i * (h_0 + yc - l_0i - u) - F_ic; % LFA-1 forces
    F2 = k_T * (h_0 - dh - l_0T - u) - F_Tn; % TCR forces

    % Return the function result
    F = [F1; F2];
end
