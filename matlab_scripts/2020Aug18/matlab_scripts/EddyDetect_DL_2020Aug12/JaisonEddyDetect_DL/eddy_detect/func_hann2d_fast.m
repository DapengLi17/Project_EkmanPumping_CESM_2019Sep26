function hsm = func_hann2d_fast ( var )
%
% function hsm = func_hann2d_fast ( var )
%
% Description : Hanning smoother/fitler for 2D-Matrix, with no
%                 data loss due to NaNs and preserves integral sum.
%                 Integral sum is preserved via treating NaN points,
%                 edges and corners accurately. Faster than func_HANN2D.m
%
%  Input  : var : 2D-matrix
%  Output : hsm : 2D-matrix, smoothed version of var 
%
%  Comments : This is a complete rewrite of original func_hann2d.m
%                  and a modified version of func_HANN2D.m
%
%  Written By : Jaison Kurian (jaison@atmos.ucla.edu)
%  Written On : Sep/14/2010
%  Reference  : Roms Tools, func_hann2d.m, func_HANN2D.m

%
%--------------------------------------------------------------------

  [jsz isz] = size(var) ;
   
  % the treament of edges and corners in func_HANN2D can be mimicked by
  %    duplicating sides on to the W,E,S,N of the respective edges and
  %    replacing corner points with NaN's. Then smoothing can be done in 
  %    a single step. 

  var_ext   = zeros (jsz+2,isz+2)       ; % add 1-more line parallell to
  var_ext(2:jsz+1,2:isz+1) = var        ; %   each of 4-sides
  var_ext(2:jsz+1,1)       = var(:,1)   ; % duplicate W-side
  var_ext(2:jsz+1,isz+2)   = var(:,isz) ; % duplicate E-side
  var_ext(1,2:isz+1)       = var(1,:)   ; % duplicate N-side
  var_ext(jsz+2,2:isz+1)   = var(jsz,:) ; % duplicate S-side
  var_ext(1,1)             = NaN        ; % NW-corner
  var_ext(1,isz+2)         = NaN        ; % NE-corner
  var_ext(jsz+2,1)         = NaN        ; % SW-corner
  var_ext(jsz+2,isz+2)     = NaN        ; % SE-corner

  % npts is used to count number of valid neighbors    
  npts               = var_ext * 0 + 1 ;
  npts(isnan(npts))  = 0 ;

  % replace NaNs with 0 to find a no-NaN sum
  var_ext(isnan(var_ext))    = 0 ;

  % initialize count and sum variables
  cc      = zeros(size(var)) ;
  varS    = zeros(size(var)) ;

  % interior points  (N,S,W,E)
  nj = jsz+2 ;
  ni = isz+2 ;
  cc =  npts(2:nj-1,2:ni-1) .* (npts(1:nj-2,2:ni-1) + ....
                                   npts(3:nj,2:ni-1) + ....
                                 npts(2:nj-1,1:ni-2) + npts(2:nj-1,3:ni) ) ;
  varS =  (var_ext(1:nj-2,2:ni-1) + ....
                                   var_ext(3:nj,2:ni-1) + ....
                                 var_ext(2:nj-1,1:ni-2) + var_ext(2:nj-1,3:ni) ) ;

  cc(cc==0) = NaN ; % bring back NaN points in original data.
  weight = 8-cc ;   % this is the weight for values on each grid point,
                    %    based on number of valid neighbours 
  hsm = ( varS + weight .* var_ext(2:jsz+1,2:isz+1)) ./8 ; % final smoothed version of var
