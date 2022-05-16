%% ========================  readme  =============================
% 
% DESCRIPTION:
% 
%  A function to compute gradient for uneven grids. 
%
% update history:
% v1.0 DL 2020Jan15
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
%   x_2d       - 2D x cartesian coordinates in E-W direction, [m]
%               e.g. x = [1 2 10;1 2 10;1 2 10]
%   y_2d       - 2D y cartesian coordinates in N-S direction, [m]
%               e.g. y = [6 6 6;4 4 4;0 0 0]
%   f          - 2D input field, e.g. sea surface height
%
% OUTPUT:
%   dfdx       - geostrophic velocity in E-W (x) direction [m/s] 
%   dfdx       - geostrophic velocity in N-S (y) direction [m/s] 
%
% EXAMPLE:
%   [dfdx,dfdy] = CalcGradient4UnevenGridsFunc(x_2d,y_2d,f)
% 
% EXTRA NOTES:
%   this function can also compute gradient for even grids if x_2d 
%   and y_2d are even spaced.
%  
% REFERENCE:
%   Matlab built-in function gradient.m  
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



