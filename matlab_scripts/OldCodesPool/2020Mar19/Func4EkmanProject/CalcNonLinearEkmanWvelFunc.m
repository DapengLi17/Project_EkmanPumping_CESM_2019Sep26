%% ========================  readme  =============================
% 
% DESCRIPTION:
% 
%  A function to compute Non-Linear Ekman vertical velocity W.
%
% update history:
% v1.0 DL 2019Sep28
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
%   rho_w     - seawater density (size: 1*1) computed from sw_dens, [kg/m3] 
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
%   W         - vertical Ekman pumping velocity [m/s] 
%
% EXTRA NOTES:
%   If use zero for the relative vorticity (kesai), the linear Ekman W vel is computed. 
%   If use nonzero for the relative vorticity, the non-liner Ekman W vel is computed 
%   using Stern (1965) equation. 
% 
%
% REFERENCE:
%   Gaube et al. 2015 JPO Satellite Observations of Mesoscale Eddy-Induced
%   Ekman Pumping, Eq (1)
%   Stern 1965 Interaction of a uniform wind stress with a geostrophic
%   vortex, Eq (1) 
% ====================================================================

function [W] = CalcNonLinearEkmanWvelFunc(rho_w,x,y,taux,tauy,f,kesai)

%% === data analysis ===
[curlz,~] = curl(x,y,taux./(f+kesai),tauy./(f+kesai)); 
W = curlz./rho_w;
W(abs(kesai./f)>0.9) = nan; % Stern Eq only apply for small Ro num 
% ======================

end
