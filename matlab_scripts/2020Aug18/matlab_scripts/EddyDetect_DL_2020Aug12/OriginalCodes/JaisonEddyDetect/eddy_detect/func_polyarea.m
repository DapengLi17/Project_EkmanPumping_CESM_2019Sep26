function parea = func_polyarea (x, y) ;
%
%   function parea = func_polyarea (x, y) ;
%
%   Description : Find the area of closed polygon, given by vertices X and Y.
%                     This is intented for the use with closed contours, where
%                     the closing is enforced. Check it manually if this function
%                     is used for any other purposes.
%
%   Input :  1. x : float, vector, X-vertices of polygon
%            2. y : float, vector, y-vertices of polygon  
%   
%   Output   1. parea ; float, scalar, area of closed polygon
%
%   
% 

    if ( nargin < 2 )
        error ('\n\t FATAL (%s) : Need two input args, X nd Y vertices of polygon.\n',mfilename)
    end
    if ( ~isvector(x) || ~isvector(y) )
        error ('\n\t FATAL (%s) : Input args (X and Y) should be vectors.\n',mfilename)
    end 

% Force closing of polygons
    
    x(end+1) = x(1) ;  % force closing of polygon
    y(end+1) = y(1) ;  %   "

% Find area of polygon

    npts     = length (x) ;
    tmp      =  ( x(1:npts-1) .* y(2:npts) ) - ( x(2:npts) .* y(1:npts-1) ) ;
    parea    = 1/2 * sum ( tmp ) ;
    parea    = abs ( parea ) ;    

% DONE

    return
