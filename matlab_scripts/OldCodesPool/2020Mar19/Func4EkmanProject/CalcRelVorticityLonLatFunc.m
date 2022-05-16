%% ========================  readme  =============================
% 
% DESCRIPTION:
% 
%  A function to compute relative vorticity kesai.
%
% update history:
% v1.0 DL 2019Oct03
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
%   u      - wind stress matrix (size: m*n) in E-W direction [N/m2]
%               e.g. taux = [0.3 0.3 0.3;0.2 0.4 0.2;0.1 0.1 0.1]
%   v      - wind stress matrix (size: m*n) in N-S direction [N/m2] 
%               e.g. tauy = -[0.1 0.2 0.3;0.4 1.0 0.6;0.7 0.8 0.9]
%
% OUTPUT:
%   kesai         - vertical Ekman pumping velocity [m/s] 
%
% EXTRA NOTES:
%   If use zero for the relative vorticity (kesai), the linear Ekman W vel is computed. 
%   If use nonzero for the relative vorticity, the non-liner Ekman W vel is computed 
%   using Stern (1965) equation. 
% 
%   This function passed test, see testCalcRelVorticityFunc.m for details
%
% REFERENCE:
% ====================================================================

function [kesai] = CalcRelVorticityLonLatFunc(lon_1d,lat_1d,u,v)

%% === data analysis ===
R       = 6371000; % [m]
dlatdy  = 1./R; % deg lat per km in y direction, dy = R.*dlat
dlondx  = 1./(R.*cosd(lat_1d)); % deg lon per km in x direction, dx = R.* cos(lat).*dlon

[~, py] = gradient(u, lon_1d, lat_1d);
[qx, ~] = gradient(v, lon_1d, lat_1d);

pvpx    = qx.*dlondx;
pupy    = py.*dlatdy;
kesai   = pvpx - pupy;
% ======================

end
