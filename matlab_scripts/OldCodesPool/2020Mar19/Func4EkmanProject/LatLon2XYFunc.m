%% ========================  readme  =============================
% 
% DESCRIPTION:
% 
%  A function to convert lat and lon to x and y.
%
% update history:
% v1.0 DL 2019Oct07
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
%   Lat_vec  - column latitude vector [deg], e.g. Lat_vec = [45 :-0.1: 27]';
%              Lat_vec(end) is the bottom point!  
%
%   Lon_vec  - row longtitude vector [deg], e.g. Lon_vec = [131:0.1:169];
%              Lon_vec(1) is the left point! 
% 
% OUTPUT:
%   x_vec    - x coordinates (vector) relative to the left bottom corner [m] 
% 
%   y_vec    - y coordinates (vector) relative to the left bottom corner [m] 
%
%   x        - meshgrid x matrix (see codes below)
%
%   y        - meshgrid y matrix (see codes below)
%
% EXTRA NOTES:
%   It calls lldistkm.m to compute the distance in km. 
%
% REFERENCE:
% ====================================================================

function [x_vec,y_vec,x,y] = LatLon2XYFunc(Lon_vec,Lat_vec)

%% === data analysis ===
% original point (left bottom corner)
latlon1 = [Lat_vec(end) Lon_vec(1)];

% compute deltaX
for i = 1 : length(Lon_vec)
    latlon2 = [Lat_vec(end) Lon_vec(i)];
    [d1x(i) ~]= lldistkm(latlon1,latlon2);
    clear latlon2
end
% diff(d1x)
x_vec = d1x.*1000; % km to m 

% compute deltaY
for i = 1 : length(Lat_vec)
    latlon2 = [Lat_vec(i) Lon_vec(1)];
    [d1y(i) ~]= lldistkm(latlon1,latlon2);
    clear latlon2
end
% diff(d1y)
y_vec = [d1y.*1000]'; % km to m 

[x,y]=meshgrid(x_vec,y_vec);
% ======================

end
