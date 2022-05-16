%% ========================  readme  =============================
% 
% DESCRIPTION:
% 
%  A function to compute z component of curl. 
%
% update history:
% v1.0 DL 2020Jan15
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
%   x_2d      - 2D x cartesian coordinates in E-W direction, [m]
%               e.g. x = [1 2 10;1 2 10;1 2 10]
%   y_2d      - 2D y cartesian coordinates in N-S direction, [m]
%               e.g. y = [6 6 6;4 4 4;0 0 0]
%   u_2d      - 2D velocity in E-W direction, [m]
%   v_2d      - 2D velocity in N-S direction, [m]
%
% OUTPUT:
%   Curlz     - z component of curl
%
% EXAMPLE:
%   Curlz = CalcCurlz4UnevenGridsFunc(x_2d,y_2d,u_2d,v_2d)
% 
% EXTRA NOTES:
%   curlz = partial v over partial x - partial u over partial y
% 
% REFERENCE:
%   N/A
% ====================================================================

function Curlz = CalcCurlz4UnevenGridsFunc(x_2d,y_2d,u_2d,v_2d)

%% === data analysis ===
  [dvdx, ~] = CalcGradient4UnevenGridsFunc(x_2d,y_2d,v_2d);
  [~, dudy] = CalcGradient4UnevenGridsFunc(x_2d,y_2d,u_2d);
  Curlz     = dvdx - dudy;  
% --------------------
% ======================

end





