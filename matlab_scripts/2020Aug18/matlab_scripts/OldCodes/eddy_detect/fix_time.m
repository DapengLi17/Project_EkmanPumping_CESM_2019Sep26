
% somehow the vale for the time axis is missing at last point. Fix it.
clear all ; close all

 etrack = 'etrack_usw51_sla_tmp.nc';

 nc = netcdf ( etrack, 'write') ;
     time = nc{'time'}(:) ;
     tlen = length(time)  ;
     dt   = time(end-1)-time(end-2) ;
     nc{'time'}(tlen) = time(end-1) + dt ;
 nc = close (nc) ;
