%%

function [x_2d,y_2d] = ConvertLonLat2XY4UnevenGridsFunc(lon_1d,lat_1d)
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

end

