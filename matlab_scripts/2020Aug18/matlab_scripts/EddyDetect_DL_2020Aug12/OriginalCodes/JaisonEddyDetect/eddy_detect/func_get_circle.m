function  [cx  cy] = func_get_circle (x0, y0, r, npts)
%
% Description : Get points on a circle, with specified (x,y) center and radius
%                    (and optional number of points too!).
%  
% Input     : 1  - x0, scalar, center X of circle
%             2  - y0, scalar, center Y of circle
%             3  - r,  scalar, radius
%             4  - npts, scalar, number of points (optional)
%
% Output    : 1  - cx, circle x-points
%             2  - cy, circle y-points
%
% Example   :  [cx cy] = func_get_circle (5, 5, 3, 256)
%              plot (cx, cy, '-')
%
% Written By : UCLA ROMS Team (jaison@atmos.ucla.edu) 
% Written On : June/05/2008
% Tool       : Eddy Tracker
%
%-----------------------------------------------------

%-----------
% Initialize   
%-----------
 
   if ( nargin < 3 ) 
       error ([mfilename ':argchk'], '\n\t FATAL (%s) : Need at least 3 inputs, x0, y0 & radius.\n',mfilename)
   elseif ( nargin == 3 )
       npts = 256;     
   end 

   if ( ~isscalar(x0) || ~isscalar(y0) || ~isscalar(r) || ~isscalar(npts))
       error ([mfilename ':argchk'], '\n\t FATAL (%s) : Need SCALAR inputs.\n',mfilename)
   end

%----------------------
% Find points on circle   
%----------------------

   
   theta = (0:npts)*2*(4*atan(1))/npts ;
   cx  = x0 + r * cos(theta) ;
   cy  = y0 + r * sin(theta) ;

%-----------
% clean up   
%-----------

   clear x0 y0 r npts theta
