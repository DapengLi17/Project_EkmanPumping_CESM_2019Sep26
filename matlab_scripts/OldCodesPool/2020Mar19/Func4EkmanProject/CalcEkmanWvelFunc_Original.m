%% ========================  readme  =============================
% 
% DESCRIPTION:
% 
%  A function to compute Ekman vertical velocity W.
%
% update history:
% v1.0 DL 2019Sep28
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
%   rho0      - seawater density (size: 1*1) computed from sw_dens, [kg/m3] 
%               e.g. rho0 = 1020; 
%   x         - x-y cartesian coordinates in E-W direction (size: m*n) [m]
%               e.g. x = [1 2 3;1 2 3;1 2 3]
%   y         - x-y cartesian coordinates in N-S direction (size: m*n) [m]
%               e.g. y = [6 6 6;4 4 4;2 2 2]
%   taux      - wind stress matrix (size: m*n) in E-W direction [N/m2]
%               e.g. taux = [0.3 0.3 0.3;0.2 0.4 0.2;0.1 0.1 0.1]
%   tauy      - wind stress matrix (size: m*n) in N-S direction [N/m2] 
%               e.g. tauy = -[0.1 0.2 0.3;0.4 1.0 0.6;0.7 0.8 0.9]
%   f         - Coriolis parameters (size: m*n) computed from sw_f, [s-1]
%               e.g. f = [3 3 3;2 2 2;1 1 1]
%   kesai     - relative vorticity (size: m*n), it can be computed from 
%               [curlz,cav] = curl(x,y,u,v); [1/s]
%               e.g. kesai = [0.9 0.8 0.7;0.6 0.5 0.4;0.3 0.2 0.1]
%
% OUTPUT:
%   W_LI      - linear Ekman pumping vertical velocity [m/s] 
%   W_NL      - nonlinear Ekman pumping vertical velocity [m/s]
%
% EXTRA NOTES: 
% 
%   This function passed test, see testCalcEkmanWvelFunc.m for details
%
% REFERENCE:
%   McGillicuddy, D. J., et al. (2008). Response to comment on 
%   "eddy/wind interactions stimulate extraordinary mid-ocean 
%   plankton blooms". science, 320(5875), 448-448.
%   Gaube et al. 2015 JPO Satellite Observations of Mesoscale Eddy-Induced
%   Ekman Pumping, paragraph around Eq (10)  
% ====================================================================

function [W_LI,W_NL] = CalcEkmanWvelFunc(rho0,x,y,taux,tauy,f,kesai)

%% === data analysis ===
[curl_tau,~]= curl(x,y,taux,tauy); 
W_LI = curl_tau./(rho0.*(f+kesai));

hx = x(1,:); 
hy = y(:,1); 
[~,PkesaiPy] = gradient(kesai, hx, hy);
[PkesaiPx,~] = gradient(kesai, hx, hy);

W_NL = (taux.*PkesaiPy - tauy.*PkesaiPx)./(rho0.*(f+kesai).^2);
% ======================

end
