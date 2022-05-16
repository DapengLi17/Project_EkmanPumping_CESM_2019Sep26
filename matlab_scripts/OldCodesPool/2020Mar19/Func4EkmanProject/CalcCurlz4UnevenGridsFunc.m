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
%   [dfdx,dfdy] = CalcGradient4UnevenGridsFunc(x_2d,y_2d,f)
% 
% EXTRA NOTES:
%   
% 
% REFERENCE:
%   
% ====================================================================

function Curlz = CalcCurlz4UnevenGridsFunc(x_2d,y_2d,u_2d,v_2d)

%% === data analysis ===
[dvdx, ~] = CalcGradient4UnevenGridsFunc(x_2d,y_2d,v_2d);
[~, dudy] = CalcGradient4UnevenGridsFunc(x_2d,y_2d,u_2d);
Curlz     = dvdx - dudy;  
% --------------------
% ======================

end




