function lon_orient = func_lon_chk(lon)
%
%    function lon_orient = func_lon_chk(lon)
%  
%  Description : Function to check the longitude orientation (e2w or w2e) of
%                  given longitude values.
%
%  Input  : lon : (float) 1D vec or 2D matrix.
%
%  Output : lon_orient :  (string) 'e2w' --> if orientation is east-to-west
%                                  'w2e' --> if orientation is west-to-east
%
%  Example : [lon_or] = func_lon_chk (lon) ;
%

%_______________________________________________________________________
%  Written By     : Jaison Kurian (jaisonk@tamu.edu
%                                  jaisonkuriann@gmail.com) 
%                   TAMU
%  Date           : Dec/01/2011
%  Reference      : func_lon_e2w.m
%  Copyright      : All Reserved
%  Modifications  : None 
%_______________________________________________________________________

%-------------
% Initialize
%-------------

    if ( nargin < 1 )
         error ( '\n\t FATAL (%s) : Need 1 Arg : lon values.\n', mfilename)
    end

    if ( ~isvector(lon))
         tmp      = lon(1,:)        ; % take all values along first lat line
         tmp(1)   = min(lon(:,1))   ; % take min of values along first lon line
         tmp(end) = max(lon(:,end)) ; % take max of values along last  lon line
         lon = tmp ;
    end
    lon_max = max(lon) ;

   if ( lon(end) < lon(1) )
         error ( '\n\t FATAL (%s) : Lon values are not monotonically increasing.\n\t\t lon_west = %g   lon_east = %g\n', mfilename, lon(1), lon(end))
    end

%%++++++++++++++++++++++++++++++++++++++++++++++
%% Check lons and find eastern and western edges
%%++++++++++++++++++++++++++++++++++++++++++++++

    % -360  -270   -180   -90     0      90    180       W
    %  +++++++++++++++++++++++++++++++++++++++++++++++++++++
    %  0     90     180   270    360    450    540       E

    % make all lon values between 0 and 360
    lon(lon<0)   = lon(lon<0) + 360.0   ; 
    lon(lon>360) = lon(lon>360) - 360.0 ; 
    if ( min(lon) < 0 || max(lon) > 360 )
        error ('\n\t FATAL (%s) : Wrong lon values : min = %g  max = %g\n',mfilename, min(lon),max(lon))
    end

    % hemisphere coverage
    hemi = 2 ; 
    if ( numel(find(lon>180)) == 0 || numel(find(lon<180)) == 0 )
       hemi = 1 ;
    end
    side = 0 ;
    if ( lon(1) < 180 && lon(end) <= 180 )  ; % lon start and end in the
       side = 1;                              %   same hemisphere
    elseif ( lon(1) >= 180 && lon(end) > 180 )%        " 
       side = 1 ;
    end 

    % orientation

    lon_orient = [] ;
    if ( lon(1) < 180 || lon(end) > 180)
        lon_orient = 'e2w' ;
    elseif (lon(1) > 180 || lon(end) < 180)
        lon_orient = 'w2e' ;
    end  
    if ( (lon(1) >= 180 && lon(end) <= 360) && lon_max < 0) % WEstern Hemi lons
        lon_orient = 'w2e' ;                                % with -ve values 
    end   
    if ( hemi == 2 && side == 1 )
        lon_orient = 'cir' ; % short for CIRCLE
    end

%-------
% DONE
%-------
    return
