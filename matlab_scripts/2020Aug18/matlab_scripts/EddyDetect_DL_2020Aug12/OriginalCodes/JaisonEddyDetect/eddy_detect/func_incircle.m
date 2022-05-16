function [in on out] = func_incircle ( x2d, y2d, x0, y0, r )
%
%    Description  :  Find indices of 2D-Matrix which lies inside, on, and/or 
%                       outside a given circle.
%
%    Input   : 1. x2d (2D-Matrix) - X-values
%              2. y2d (2D-Matrix) - Y-values
%              3. x0  (scalar)    - Ceter X of circle
%              4. y0  (scalar)    - Ceter Y of circle
%              5. r   (scalar)    - Radius of circle
%
%    Output  : 1. in  (2D-Matrix) - 1 if inside  the circle, else 0
%              2. on  (2D-Matrix) - 1 if on      the circle, else 0
%              3. out (2D-Matrix) - 1 if outside the circle, else 0
%
%
%    Written By : UCLA ROMS Team (jaison@atmos.ucla.edu)
%    Written On : June/12/2008 
%    Copyright @ UCLA ROMS Team
%    Tool       : Eddy Tracker
%    Version    : 1.0
%
%------------------------------------------------------------------------

     dist = sqrt ( (x2d - x0 ).^2 + (y2d - y0 ).^2 ) ;
     
     in   = dist < r ;
     on   = abs(dist - r) < 1e-12 ;
     out  = dist > r ;

     return
