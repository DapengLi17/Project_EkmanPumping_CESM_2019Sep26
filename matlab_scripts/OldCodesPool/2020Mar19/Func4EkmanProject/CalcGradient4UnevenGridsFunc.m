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

function [dfdx,dfdy] = CalcGradient4UnevenGridsFunc(x_2d,y_2d,f)

%% === data analysis ===
[m,n] = size(f);

% --- compute dfdx ---
% Take forward differences on left and right edges
dfdx(:,1) = (f(:,2) - f(:,1))   ./(x_2d(:,2)-x_2d(:,1));
dfdx(:,n) = (f(:,n) - f(:,n-1)) ./(x_2d(:,n)-x_2d(:,n-1));

% Take centered differences on interior points
dfdx(:,2:n-1) = (f(:,3:n) - f(:,1:n-2)) ./ (x_2d(:,3:n) - x_2d(:,1:n-2));
% --------------------

% --- compute dfdy ---
% Take forward differences on top and bottom edges
dfdy(1,:) = (f(2,:) - f(1,:)) ./ (y_2d(2,:) - y_2d(1,:));
dfdy(m,:) = (f(m,:) - f(m-1,:)) ./ (y_2d(m,:) - y_2d(m-1,:));

% Take centered differences on interior points
dfdy(2:m-1,:) = (f(3:m,:) - f(1:m-2,:)) ./ (y_2d(3:m,:) - y_2d(1:m-2,:));
% --------------------
% ======================

end



