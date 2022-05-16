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
%   ssh       - sea surface height, [m]
%   x_2d      - 2D x cartesian coordinates in E-W direction, [m]
%               e.g. x = [1 2 10;1 2 10;1 2 10]
%   y_2d      - 2D y cartesian coordinates in N-S direction, [m]
%               e.g. y = [6 6 6;4 4 4;0 0 0]
%   lat_1d    - 1D latitude array, [deg]
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

function [u,v] = CalcGeostrophyVel4UnevenGridsFunc(ssh,x_2d,y_2d,lat_1d)

%% === data analysis ===
  g = 9.8; % [m/s^2]
  f = sw_f(lat_1d);
  f_rp = repmat(f,1,size(x_2d,2)); % rp: repmat  
  
  [dsshdx,dsshdy] = CalcGradient4UnevenGridsFunc(x_2d,y_2d,ssh);
   
  u = -g./f_rp.*dsshdy;
  v = g./f_rp.*dsshdx;
% ======================

end
