clear all
clc
%
% LFA-1-ICAM-1 Catch Bond Curve Fitting using least squares method
%  
% using two-pathway model 
% Reference papers: C.E. Chan and D.J. Odde, Science, 2008; 
%                   W.E. Thomas et al, Annual Review Biophysics, 2008;
%                   V. Luca et al, Science 2017;
%                   J. Brockman et al, Frontiers in Physics 2019
%                   W. Chen 4t al, Journal of Biological Chemistry 2010;
%
% Inputs:
%
% Two-pathway model: k_off(F) = k_c0*exp(-F/F_c) + k_s0*exp(-F/F_s)
%
% Single molecule 2D micripipette LFA-1-ICAM-1 data from V. Luca et al Fig. 1A, Ca2+/Mg2+/CXCL12;
% 
% f_LFA-1-ICAM-1 = [0, 4, 7.5, 12.5, 18, 24] pN
% tau_LFA-1-ICAM-1 = 1/k_off(F) = [0.4, 0.75, 1.125, 1.5, 0.6, 0.25] s
%
% f_TCR-n-pMHC = [0, 4, 7.5, 12.5, 18, 24] pN
% tau_TCR-n-pMHC = 1/k_off(F) = [0.4, 0.75, 1.125, 1.5, 0.6, 0.25] s
% f_TCR-c-pMHC = []
% f_BOND is the applied load (pN)
% tau_BOND is the dwell time in the bound state (s) w.r.t. f_BOND
% Outputs:
%
% k_c0, k_s0, F_c, F_s all determined a priori from model fit
% A = [k_c0, k_s0, F_c, F_s] <- parameter extraction passed into vector A
%
% data points to fit from literature
f_LFA1ICAM1 = [0, 4, 7.5, 12.5, 18, 24]; % pN
tau_LFA1ICAM1 = [0.4, 0.75, 1.125, 1.5, 0.6, 0.25]; % s

% Fit model using least squares method
tau_fitICAM1 = @(A, f_LFA1ICAM1) 1./(A(1).*exp(-f_LFA1ICAM1./A(2))+...
                                A(3).*exp(f_LFA1ICAM1./A(4)));
A0 = [10, 4, 1, 2]; % arbitrarily determined initial conditions for fitting                            
A = lsqcurvefit(tau_fitICAM1, A0, f_LFA1ICAM1, tau_LFA1ICAM1);
f = linspace(f_LFA1ICAM1(1),f_LFA1ICAM1(end));

figure(1)
plot(f_LFA1ICAM1,tau_LFA1ICAM1,'ko',f,tau_fitICAM1(A,f),'b-');
title('Bond Lifetime vs. Force LFA-1-ICAM-1');
xlabel('Force (pN)'); ylabel('Bond Lifetime (s)');