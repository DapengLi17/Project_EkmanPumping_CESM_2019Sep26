%% ========================  readme  =============================
% 
% DESCRIPTION:
% 
%  A function to compute geostrophic velocity based on sea surface height. 
%
% update history:
% v1.0 DL 2020Jan15
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
%   x         - x-y cartesian coordinates in E-W direction (size: m*n) [m]
%               e.g. x = [1 2 3;1 2 3;1 2 3]
%   y         - x-y cartesian coordinates in N-S direction (size: m*n) [m]
%               e.g. y = [6 6 6;4 4 4;2 2 2]
%   eta       - sea surface height seawater density (size: 1*1) computed from sw_dens, [m] 
%               e.g. rho0 = 1020; 
%   f         - Coriolis parameters (size: m*n) computed from sw_f, [s-1]
%               e.g. f = [3 3 3;2 2 2;1 1 1]
%   g         - acceleration due to gravity, [m/s^2]
%
% OUTPUT:
%   u         - geostrophic velocity in E-W (x) direction [m/s] 
%   v         - geostrophic velocity in N-S (y) direction [m/s] 
%
% EXAMPLE:
%   [u,v] = CalcGeostrophicVelFunc(x,y,eta,f,g);
% 
% EXTRA NOTES:
%   Geostrophic velocity (u,v) are formulated as:
%   u = -g/f partial(eta) over partial(y)
%   v = g/f partial(eta) over partial(x)
% 
% REFERENCE:
%   Siegelman et al. 2020 Nature(geoscience) Enhanced upward heat transport
%   at deep submesoscale ocean fronts, Eq (1) in Methods
%   Park 2004 JGR Determination of the surface geostrophic velocity field
%   from satellite altimetry, Eq(2) in Method
% ====================================================================

function [u,v] = CalcGeostrophVelLonLatFunc(lon_1d,lat_1d,eta)

%% === data analysis ===
g = 9.8; % [m/s^2]
R = 6371000; % [m]
f = sw_f(lat_1d);
f_rp = repmat(f,1,length(lon_1d)); % rp: repmat

dlatdy = 1./R; % deg lat per km in y direction, dy = R.*dlat
dlondx = 1./(R.*cosd(lat_1d)); % deg lon per km in x direction, dx = R.* cos(lat).*dlon
dlondx_rp = repmat(dlondx,1,length(lon_1d));     

    [px, py] = gradient(eta, lon_1d, lat_1d);
   
    u = -g./f.*py.*dlatdy;
    v = g./f.*px.*dlondx_rp;
% ======================

end
