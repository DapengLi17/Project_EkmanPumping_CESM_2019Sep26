function [uniq_eddy dist indx] = func_uniq_eddy (cx, cy, cr, peddy)
%
%  Description : Matlab function to determine whether an eddy for a specific
%                  contour value with its center at cx,cy has been already
%                  identified (in edx,edy,edb) or not (circle fit method).
%
%  Input : cx    - (float, scalar, km)  center x-position of eddy to be tested
%          cy    - (float, scalar, km)  center y-position of eddy to be tested
%          cr    - (float, scalar, km)  radius of eddy to be tested    
%          peddy - (struct, float)      cyc or acyc structure from func_get_eddies
%                              which contains the following fields:
%                           peddy.xkm
%                           peddy.ykm
%                           peddy.rkm
% 
%  Output : uniq_eddy (float, scalar),  1 = eddy is uniq
%                                       0 = eddy has been already identified 
%           dist                     , distance of current eddy center from nearest
%                                            eddy center in previous snapshot
%           indx                     , grid (j,i) index correspoinding to above dist
%
%  Method : - find the distance (in km) between the center position of current
%                    eddy and all previously identified eddies. 
%           - if the distance is </= maximum value of radius between eddis being
%                    compared, then the current eddy is not a unique one.

%  Written By : Jaison Kurian (jaison@atmos.ucla.edu)
%  Written On : Sep/22/2008 (for ellipse case)
%  Modifications :  Oct/28/2008 Adapted for circle fitting case
%                   Dec/02/2008 Distace vs radius test changed to suit practical
%                                 situations.
% 
%  Part of Eddy tracking toolkit
%
%---------------------------------------------------------------------------

%-------------        
% Initialize
%-------------        

    % argument checks

    if ( nargin  < 4 ) 
       error([mfilename ':argchk'], '\n   FATAL (%s) : Too few input arguments. \n', mfilename)
    end

    if ( ~isscalar (cx) || ~isscalar (cy) || ~isscalar(cr) ||~isstruct (peddy) ) 
       error([mfilename ':argchk'], '\n   FATAL (%s) : Wrong input variable type: Need scalar, scalar, scalar & struct. \n', mfilename)
    end

%----------------------------
% Check if the eddy is unique
%----------------------------

    % distance and index of the closest eddy
    [dist indx]  = func_NNR_indx_dist ( peddy.xkm,   peddy.ykm,   cx,  cy ) ;

    % circle parameters of the closest eddy
    px = peddy.xkm(indx) ; % x-position (in km) 
    py = peddy.ykm(indx) ; % y-position (in km) 
    pr = peddy.rkm(indx) ; % radius (in km) 
    
    % distance between current eddy and closest eddy in previous iteration
    dist = sqrt ( (cx - px).^2 + (cy - py).^2 ) ;

    % uniq eddy
    uniq_eddy = 0             ; % assume it is not a unique eddy
    if ( dist > cr+pr+0.05*min(cr, pr) ) ; % if centers are separated by 
                                           %  sum_of(radii + 5%_of_small_radius)
         uniq_eddy = 1        ; %    then it is a unique eddy
    end 
     
%-----
% Done
%-----

    clear pr px py 
    return
