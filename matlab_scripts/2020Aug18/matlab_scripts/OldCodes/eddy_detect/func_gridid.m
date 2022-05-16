function ids = func_gridid ( xgrid, ygrid, x, y )
%
% Description : Find the [row, colum] index in a given grid (specified by xgrid/ygrid), closest to the given x and y values. Grids should have unique
%                  monotonous values, without repetition.
%
% Input : 1. xgrid (2D-Matrix)  - grid X values, all positive, monotonic
%         2. ygrid (2D-Matrix)  - grid Y values, all positive, monotonic
%         3. xval (scalar)      - X value of point
%         4. yval (scalar)      - Y value of point
%
% Output : 1. ids (vector) - [row column] 2-elemant vector for the grid point
%                                close to (x,y)
%
% Written By : Jaison Kurian (jaison@atmos.ucla.edu)
% Written On : June/23/2008
%
%--------------------------------------------------------------------------

   if ( nargin < 4 || ~func_shape(xgrid,2) || ~func_shape(ygrid,2) || ~func_shape(x,0) || ~func_shape(y,0)) 
        error ([mfilename ':argchk'], '\n FATAL (%) : Wrong number/type of inputs: need xgrid & ygrid (2D-Matrices) and x & y values (scalars).', mfilename) 
   end 

   xdif = xgrid - x ;
   ydif = ygrid - y ;
   tmp  = abs ( xdif ) + abs ( ydif ) ;
   [tmp1 ids] = min2d ( tmp ) ;  

   return 
