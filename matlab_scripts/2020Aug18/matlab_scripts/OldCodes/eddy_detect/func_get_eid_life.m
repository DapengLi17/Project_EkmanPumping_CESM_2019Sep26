function [eid edb edd] = func_get_eid_minlife (eddi, min_life)
%
%  Description : Find idices (with respect to eddy tracking output
%                   NetCDF file) of eddies with life time => min_life
%                   (in snapshots).
%
%  Input  :  eddi - float, 2D-Matrix (time, eid) any eddy variable, properly 
%                     masked (NaNs applied) (often output from func_init_eddy)
%            min_life - float, scalar, minimum snapshots required 
%
%  Output :  eid's of selected eddies
%            b_id   - linear index of birth time
%            d_id   - linear index of death time 
%
%  Ref : func_pick_eddy.m
%
%  Jaison Kurian
%  June/23/2009
%
% ----------------------------------------------------------------------------

     if ( nargin < 2 )
        error('\n   FATAL (%s) : Wrong number of input arguments (%4i). Need 2.\n', mfilename, nargin)
     end

     if ( ~func_shape(eddi,2) )
        error('\n   FATAL (%s) : Wrong dimensions for arg1, need 2D-Matrix (time, eid).\n', mfilename)
     end
     if ( ~isscalar(min_life) )
        error('\n   FATAL (%s) : Wrong dimensions for arg2, need scalar.\n', mfilename)
     end

     [tsize isize] = size(eddi) ;
     if ( min_life >= tsize )
        error('\n   FATAL (%s) : min_life (arg2, %g) > time size (%g).\n', mfilename, min_life, tsize)
     end

     eid         =  [] ;

     ref = eddi ;
     ref( isnan(ref) ) = 0   ;            % set NaN values to zero, for reference variable
     isz = find ( sum(ref,1),1,'last' ) ; % Find the I/column-index for the last valid eddy.
     if ( isz < isize )                   % if it is less than the X/column-size 
        ref        = ref(:,1:isz) ;       % re-size all variables used below
        eddi       = eddi(:,1:isz) ;      %         "
     end

     % Find eddy date/time of eddy birth and death

     ref       = ref'     ;               % transposing is just to use an already existing 
                                          %    algorithm 
     ref_li    = zeros( [size(ref,1), 1] ) ; % init column vector
     ref_lf    = zeros( [size(ref,1), 1] ) ; % init   "
     [r c]     = find (ref) ;                % indices of non-zero values
     ref_lf(r) = c ;                         % eddy end   T-ids (column index of last  non-zero value)
     ref_li(r(end:-1:1)) = c(end:-1:1) ;     % eddy start T-ids (column index of first non-zero value)      

     ref_lf    = ref_lf' ; % transpose back to original shape
     ref_li    = ref_li' ; % transpose back to original shape 

     mask_life  = ref_lf - ref_li + 1 ;
     mask_life ( mask_life < min_life ) = 0 ;
     eid        = find (mask_life) ;

     % find indices corresponding to eddy birth/death time

     edb     = ref_li(eid) ;
     edd     = ref_lf(eid) ;

return
