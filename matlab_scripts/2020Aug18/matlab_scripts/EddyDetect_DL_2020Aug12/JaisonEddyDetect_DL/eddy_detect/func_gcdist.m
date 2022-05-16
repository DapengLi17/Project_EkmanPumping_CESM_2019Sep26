function distkm = func_gc_dist (lat1, lon1, lat2, lon2, strict)
%
% function distkm = func_gc_dist (lat1, lon1, lat2, lon2, strict)
%
% Description : Find distance between two points on a sphere using 
%                  the great-circle-distancemethod.
%
% Test case   : http://www.pmel.noaa.gov/maillists/tmap/ferret_users/fu_2005/msg00011.html
%                  Find the distance between Seattle (47.55N, 122.33W) and Tokyo
%                  (35.75N, 139.5E):   7707.6 km
%
%  Inputs : 1 - lat1 (float, any dim, degrees N) latitude/s  of point/array 1
%           2 - lon1 (float, any dim, degrees N) longitude/s of point/array 1
%           3 - lat2 (float, any dim, degrees E) latitude/s  of point/array 2
%           4 - lon2 (float, any dim, degrees E) longitude/s of point/array 2
%           5 - strict (logical) : lon1 and lon2 should be mutually inclusive
%                                  (over same region) 
%
%               For 2 and 4, Western longitudes are OK (will be converted to E).
%
%  Outputs : distance in km.

%________________________________________________________________________________
%  Written By     : Jaison Kurian (jaisonk@tamu.edu
%                                  jaisonkuriann@gmail.com)
%                   Oceanography Dept., TAMU
%  Date           : Feb/08/2012
%  eTRACK Version : 2.00
%  Ref. Version   :                    
%  Copyright      : See etrack_manual.pdf
%  Modifications  : see below
%_______________________________________________________________________
%  Feb/08/2013    :  Added lon check using func_lon_chk Feb/08/2012
%_______________________________________________________________________


%----------
% Iniialize
%----------

    if ( nargin == 4 )
       strict = 0 ;
    end
    if ( strict )
       if ( numel(lon1) == 1)
           if ( min(min(lon2)) > lon1 || max(max(lon2)) < lon1 )
              error ('\n\t FATAL (%s) : Lon1(arg2) is outside Lon2(arg4) range.\n\t\t Lon1 = %6.2f  Lon2 = %6.2f to %6.2f',mfilename,lon1,min(min(lon2)),max(max(lon2)))
           end
       end
       if ( numel(lon2) == 1)
           if ( min(min(lon1)) > lon2 || max(max(lon1)) < lon2 )
              error ('\n\t FATAL (%s) : Lon2(arg4) is outside Lon1(arg2) range.\n',mfilename)
           end
       end
    end

    % convert Western Longitudes to Eastern

    if ( numel(lon1) > 1 )
        or1 = func_lon_chk(lon1) ;
        if ( strcmp(or1,'e2w') == 1 )
           lon1(lon1 < 0) = lon1(lon1 < 0) + 360 ;
        end 
    end
    if ( numel(lon2) > 1 )
        or2 = func_lon_chk(lon2) ;
        if ( strcmp(or2,'e2w') == 1 )
           lon2(lon2 < 0) = lon2(lon2 < 0) + 360 ;
        end 
    end

    % constants

    pii       = 4 * atan (1)  ;
    deg2rad   = pii/180       ; % radians 
    rad_Earth = 6371.2        ; % km, mean radius of earth
    circumf   = 2 * pii * rad_Earth ;
 
%----- 
% Main
%----- 

    a = sin((lat2-lat1)*deg2rad/2).^2 + cos(lat1*deg2rad) .* cos(lat2*deg2rad) .* sin((lon2-lon1)*deg2rad/2).^2;
    c = 2 * atan2(sqrt(a),sqrt(1 - a)) ;
    distkm = rad_Earth * c ;

    return
