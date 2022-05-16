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
%   W_TE      - total Ekman pumping W vel from Stern (1950) [m/s]
%   W_LE      - linear Ekman pumping vertical velocity [m/s] 
%   W_NE      - nonlinear Ekman pumping vertical velocity [m/s]
%
% EXTRA NOTES: 
% 
%   This function passed test by DL on 2020Mar17, see 
%   testCalcEkmanWvel4UnevenGridsFunc.m for details.
%
% REFERENCE:
%   McGillicuddy, D. J., et al. (2008). Response to comment on 
%   "eddy/wind interactions stimulate extraordinary mid-ocean 
%   plankton blooms". science, 320(5875), 448-448.
%   Gaube et al. 2015 JPO Satellite Observations of Mesoscale Eddy-Induced
%   Ekman Pumping, paragraph around Eq (10)  
% ====================================================================

function [W_TE,W_LE,W_NE] = CalcEkmanWvel4UnevenGridsFunc(rho_w,x_2d,y_2d,taux,tauy,f_2d,kesai)

%% === data analysis ===

  VorTotal = f_2d+kesai;

  Curlz4TE = CalcCurlz4UnevenGridsFunc(x_2d,y_2d,taux./VorTotal,tauy./VorTotal); 
  W_TE = Curlz4TE./rho_w;

  Curlz4LE = CalcCurlz4UnevenGridsFunc(x_2d,y_2d,taux,tauy); 
  W_LE = Curlz4LE./VorTotal./rho_w;

  [dVorTotaldx, dVorTotaldy] = CalcGradient4UnevenGridsFunc(x_2d,y_2d,VorTotal);
  W_NE = (taux.*dVorTotaldy - tauy.*dVorTotaldx)./(VorTotal.^2)./rho_w;

% discard Wvel if (1) near equator (f<f_min)
%  Wvel(abs(f_2d)<f_min) = nan; 
% (2) kesai ~ f (Stern Eq only apply for small Ro num)
% Wvel(abs(kesai./f_2d)>0.5) = nan;  
% ======================

end
