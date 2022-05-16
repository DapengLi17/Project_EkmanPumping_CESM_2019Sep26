
%% ========================  readme  =============================
% 
% DESCRIPTION:
% 
%  A function to compute time average for 3D array (lat,lon,time).
%
% update history:
% v1.0 DL 2019Nov14
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
%   w_raw  - 3D Wvel matrix, 1st Dimension: lat, 2nd Dimension: lon, 3rd
%               Dimension: time, [m/s]
%   T_raw  - 3D T matrix, 1st Dimension: lat, 2nd Dimension: lon, 3rd
%               Dimension: time, [deg C]
%   rho_w  - reference seawater density [kg/m3]
% 
% OUTPUT:
%   JqEddy - 3D dat matrix, vertical eddy heat flux [W/m2]
%
% EXTRA NOTES:
% N/A 
% 
% REFERENCE:
% N/A
% ====================================================================

function [JqEddy,w_prime,T_prime] = CalcVerticalEddyHeatFluxFunc(w_raw,T_raw,rho_w)

%% === data analysis ===
  Cp = 3.99e3; % [J/kg/k] the specific heat at constant pressure (Thorpe 2007) 
 
  len_time= size(w_raw,3); % length of time

  w_av    = nanmean(w_raw,3);
  w_av_3D   = repmat(w_av,1,1,len_time); % duplicate mean field to 3D
  w_prime = w_raw - w_av_3D;
  T_av    = nanmean(T_raw,3);
  T_av_3D   = repmat(T_av,1,1,len_time); % duplicate mean field to 3D
  T_prime = T_raw - T_av_3D;

  JqEddy  = nanmean(rho_w.*Cp.*w_prime.*T_prime,3);
% ======================

end
