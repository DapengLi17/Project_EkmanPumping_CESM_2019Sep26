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
%   [pFpx,pFpy,x_2d,y_2d] = CalcGradientLonLatFunc(lon_1d,lat_1d,F)
% 
% EXTRA NOTES:
%   
% 
% REFERENCE:
%   
% ====================================================================

function [pFpx,pFpy,x_2d,y_2d] = CalcGradientLonLatFunc(lon_1d,lat_1d,F)

%% === data analysis ===

R = 6371000; % Earth radius [m]

% --- convert lon and lat to x and y ---
lon_reso_rad = (lon_1d(2)-lon_1d(1))./180.*pi;
lat_reso_rad = (lat_1d(2)-lat_1d(1))./180.*pi;

dy = lat_reso_rad.*R;
y_1d = [0 : dy : (length(lat_1d)-1).*dy]'; % y=0 corresponds to the 1st element of lat_1d
y_2d = repmat(y_1d,1,length(lon_1d));

for ilat = 1 : length(lat_1d)
   dx(ilat)     = lon_reso_rad.*R.*cosd(lat_1d(ilat));
   x_2d(ilat,:) = [0 : dx(ilat) : (length(lon_1d)-1).*dx(ilat)]; % x=0 corresponds to the 1st element of lon_1d
end

% --- compute gradients ---
[dx,~]=gradient(x_2d);
[~,dy]=gradient(y_2d);
[dF_x,dF_y]=gradient(F);
dFdx=dF_x./dx;
dFdy=dF_y./dy;
% ======================

end

