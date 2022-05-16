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
%   rho_w     - seawater density (size: 1*1) computed from sw_dens, [kg/m3] 
%               e.g. rho0 = 1020; 
%   x_2d      - 2D x cartesian coordinates in E-W direction, [m]
%               e.g. x = [1 2 10;1 2 10;1 2 10]
%   y_2d      - 2D y cartesian coordinates in N-S direction, [m]
%               e.g. y = [6 6 6;4 4 4;0 0 0]
%   taux      - wind stress matrix in E-W direction [N/m2]
%               e.g. taux = [0.3 0.3 0.3;0.2 0.4 0.2;0.1 0.1 0.1]
%   tauy      - wind stress matrix in N-S direction [N/m2] 
%               e.g. tauy = -[0.1 0.2 0.3;0.4 1.0 0.6;0.7 0.8 0.9]
%   f_2d      - 2D Coriolis parameters computed from sw_f and repmat, [s-1]
%               e.g. f = [3 3 3;2 2 2;1 1 1]
%   kesai     - relative vorticity, [s-1]
%               e.g. kesai = [0.9 0.8 0.7;0.6 0.5 0.4;0.3 0.2 0.1]
%
% OUTPUT:
%   Wvel      - Ekman pumping vertical velocity [m/s] 
%   if kesai = 0, Wvel is linear Ekman vertical velocity 
%   if kesai /=0, Wvel is nonlinear Ekman vertical velocity
%
% EXTRA NOTES: 
%   Wvel is NaN near equator because f=0 at the equator
%
% REFERENCE:
%   McGillicuddy, D. J., et al. (2008). Response to comment on 
%   "eddy/wind interactions stimulate extraordinary mid-ocean 
%   plankton blooms". science, 320(5875), 448-448.
%   Gaube et al. 2015 JPO Satellite Observations of Mesoscale Eddy-Induced
%   Ekman Pumping, paragraph around Eq (10)  
% ====================================================================

function Wvel = CalcEkmanWvelFunc(rho_w,x_2d,y_2d,taux,tauy,f_2d,kesai,f_min)

%% === data analysis ===
  Curlz = CalcCurlz4UnevenGridsFunc(x_2d,y_2d,taux./(f_2d+kesai),tauy./(f_2d+kesai)); 
  Wvel = Curlz./rho_w;

% discard Wvel if (1) near equator (f<f_min)
  Wvel(abs(f_2d)<f_min) = nan; 
% (2) kesai ~ f (Stern Eq only apply for small Ro num)
% Wvel(abs(kesai./f_2d)>0.5) = nan;  
% ======================

end
