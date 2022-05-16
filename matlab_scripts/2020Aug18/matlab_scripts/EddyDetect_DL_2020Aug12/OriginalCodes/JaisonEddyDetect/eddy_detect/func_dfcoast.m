function  dist = func_dfcoast ( xrow, mrow, iid )
%
%     Description : Find the distance of given grid point (idkm) from 
%                       model's coastline.
%
%     Input : 1. xrow (1D-vector) - X-values of grid in km, along row or j=cent_j
%             2. mrow (1D-vector) - mask for row j=cent_j
%             3. iid (scalar)     - i or column value of eddy center
%
%     Output : dist (scalar) - distance in km, from coast.
%
%     Assumptions : 1. Model domain is for an eastern-boundary region, with coast
%                         to the east.
%                   2. While searching along j-rows, backward from last i-column,
%                         the grid point immediately to the east of first ocean 
%                         grid point is the coast.
%
%                                         |coast
%                                         |
%                          x    x    x   NaN   NaN   NaN    NaN
%                     (cent_j,            |                (cent_j,i=end)
%                       cent_i)
%
%    Method : Find the minimum value and index for max(xrow) - xrow_masked
%                 coast_i = above_minium_i + 1
%
%    Written By : UCLA ROMS Team (jaison@atmos.ucla.edu)
%    Written On : June/17/2008
%
%----------------------------------------------------------------------------

%-----------
% Initialize
%-----------

    if ( nargin < 3 )
        error ([mfilename ':argchk'], '\n\t FATAL (%s) : Need 2 input arguments (grid, center_ji).\n', mfilename)
    end      

%-----------------------------
% Find the distance from coast
%-----------------------------

    tmp1 = xrow .* mrow ;
    [tmp2 idm] = min ( max(xrow) - tmp1 ) ;

    ilen = numel (xrow) ; 
    idm  = min ( idm+1, ilen ) ; % safe guard

    dist = xrow(idm) - xrow(iid) ;
 
    return   
