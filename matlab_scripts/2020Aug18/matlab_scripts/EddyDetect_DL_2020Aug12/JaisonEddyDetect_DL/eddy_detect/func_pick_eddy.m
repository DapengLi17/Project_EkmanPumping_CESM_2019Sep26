function [eid estat ebx eby edx edy ebt elife] = func_pick_eddy ( eddy, min_snap, max_dfco, min_lat, max_lat, tids, tide )
%
%  Description : Matlab function to pick eddies which satisfy given criterion,
%                   from eddy tracking outputs. Criterion are following :
%
%                   maximum distance from coast : for eddy birth location                  
%                   minimum latitude            : for eddy birth location                  
%                   maximum latitude            : for eddy birth location                  
%                   time period                 : for eddy birth time
%                   minimum life/snap           : minimum life time for eddies 
%                                                   which meets all of above 4
%                                                   criterion.   
%
%  Inputs  : 1. eddy  - (mixed, struct) eddy variables from func_init_eddy
%            2. min_snap - (float, scalar) minimum lifetime (in number of 
%                              snapshots/netcdf records) required for an eddy
%            3. max_dfco - (float, scalar) maximim distance (km) from coast, 
%                                (applied to *_dcokm in eddy tracking file).
%            4. min_lat  - (float, scalar) min Lat (in Deg_N) of interested region
%            5. max_lat  - (float, scalar) max Lat (in Deg_N) of interested region
%            6. tids  - (float, scalar, optional)  starting time index to 
%                                extract data. If not specified, then set to 1.
%            7. tide  - (float, scalar, optional)  ending   time index to
%                                extract data. If not specified, then set to max.
%    
%  Outputs : 1. eid  - (float, vector) column ids of eddies which satisfies all
%                         5 criterions.
%            2. stat - (float, struct) statistics about howmany eddies has 
%                         passed each of the 5 tests.
%            3. ebx  - (float, vector) x-birth position of selected eddies (deg_E)
%            4. eby  - (float, vector) y-birth position of selected eddies (deg_N)
%            5. edx  - (float, vector) x-death position of selected eddies (deg_E)
%            6. edy  - (float, vector) y-death position of selected eddies (deg_N)
%            7. ebt  - (float, vector) birth time-index of selected eddies
%            8. elife - (float, vector) life period of selected eddies, in number
%                                of snapshots.
%

%_____________________________________________________________________________
% This is a part of Eddy Tracking Tool Kit, Copyright @ UCLA ROMS Team, 2008
% Written By : Jaison Kurian (jaison@atmos.ucla.edu)
% Written On : Nov/05/2008
%_____________________________________________________________________________

 
% Initialize

     if ( nargin < 5 || nargin == 6 )
        error([mfilename ':argchk'], '\n   FATAL (%s) : Wrong number of input arguments (%g). Need either 6 or 8.\n', mfilename, nargin)
     end  

     if ( ~isstruct(eddy) ) 
        error([mfilename ':argchk'], '\n   FATAL (%s) : Wrong input type (arg2) : should be a struct \n\t\t created b calling func_init_eddy.\n', mfilename)
     end

     [tsize isize] = size(eddy.i) ;

     if ( nargin > 5 )
        if ( ~isfloat(tids) || ~isfloat(tide) || ~func_shape(tids,0) || ~func_shape(tide,0) )
           error([mfilename ':argchk'], '\n   FATAL (%s) : Wrong arg type (args 7 & 8) : Need float scalars, \n\t\t starting and ending time indices).\n', mfilename)
        end 
        if ( min(tids,tide) < 1 || max(tids,tide) > tsize )
           error([mfilename ':argchk'], '\n   FATAL (%s) : Wrong time range (args 7 & 8) : should be \n\t\t between 1 and %g (present range %g:%g).\n', mfilename, tsize, tids, tide)
        end
     else 
        tids = 1 ;
        tide = tsize ;  
     end

     max_snap = tide-tids+1 ;
     if ( min_snap < 1 || min_snap >= max_snap )
        error([mfilename ':argchk'], '\n   FATAL (%s) : Wrong minimum snap/life (%g, arg3) \n\t\t should be 1 <= min_snap < %g.\n', mfilename, min_snap, max_snap)
        
     end

     if ( max_dfco < 20 ) 
        error([mfilename ':argchk'], '\n   FATAL (%s) : Wrong maximum distance from coast (%g, arg4) \n\t\t should be >= 20 km.\n', mfilename, max_dfco)
     end

     if ( min_lat > max_lat) 
        error([mfilename ':argchk'], '\n   FATAL (%s) : Minimum Lat (%g) > Maximum Lat (%g).\n', mfilename, min_lat, max_lat)
     end

     % init output variables with empty fields

     eid         =  [] ; 
     estat.ndco  =  [] ;
     estat.nlat  =  [] ;
     estat.nbt   =  [] ;
     estat.nlife =  [] ;
     estat.n     =  [] ;
     ebx         =  [] ;
     eby         =  [] ;
     edx         =  [] ;
     edy         =  [] ;
     ebt         =  [] ;
     elife       =  [] ;

         

%--------------------------
% Mask with given criterion
%--------------------------

     % remove the padded missing values at the end (this has been done to have a common
     %     X-dimension/axis for the eddy tracking output NetCDF file).

     clear ref ;
     ref = eddy.i ;
     ref( isnan(ref) ) = 0   ;            % set NaN values to zero, for reference variable
     isz = find ( sum(ref,1),1,'last' ) ; % Find the I/column-index for the last valid eddy.
     if ( isz < isize )                   % if it is less than the X/column-size 
        ref        = ref(:,1:isz) ;       % re-size all variables used below
        eddy.i     = eddy.i(:,1:isz) ;    %         "
        eddy.dcokm = eddy.dcokm(:,1:isz) ;%         "
        eddy.y     = eddy.y(:,1:isz) ;    %         "
        eddy.x     = eddy.x(:,1:isz) ;    %         "
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

     % make a mask for the death time, if it is the last or second last time point

     nod_byend = ref_lf ;  % no-death-by-end-of-model-run
     nod_byend (nod_byend>=tsize-1) = NaN ;
     nod_byend = nod_byend * 0 + 1 ; % make it a mask, 1D-vector, length = number of eddies

     % find indices corresponding to eddy birth time

     b_id     = sub2ind(size(eddy.i),ref_li,[1:1:size(eddy.i,2)]) ;
     d_id     = sub2ind(size(eddy.i),ref_lf,[1:1:size(eddy.i,2)]) ;
     
     % Mask 1 : Distance from the coast, at birth time.  
     %          Use NaN to mask insted of 0 as few of the eddy variables  
     %          can have a 0 value, which is valid in practical sense.

     mask_dco   = eddy.dcokm(b_id) ;
     mask_dco ( mask_dco > max_dfco ) = NaN ; 
     estat.ndco = size(mask_dco(~isnan(mask_dco))) ; % no. of eddies which has passed
                                                     %   dfco test

     % Mask 2 : Latitude

     mask_lat   = eddy.y(b_id) ;
     mask_lat (mask_lat < min_lat ) = NaN  ;
     mask_lat (mask_lat > max_lat ) = NaN  ;
     estat.nlat = size(mask_lat(~isnan(mask_lat)),2) ; 
 
     % Mask 3 : Birth time : take those eddies which have borned on/after tids (!!)
                              % and on/before tide
     mask_bt    = ref_li ;
     mask_bt ( mask_bt < tids ) = NaN ;
     mask_bt ( mask_bt > tide ) = NaN ;
     estat.nbt  = size(mask_bt(~isnan(mask_bt)),2) ; 

     % Mask 4 : Life period

     mask_life  = ref_lf - ref_li + 1 ;
     mask_life ( mask_life < min_snap ) = NaN ; 
     estat.nlife= size(mask_life(~isnan(mask_life)),2) ; 
 
     % Mask ALL

     mask       = mask_dco + mask_lat + mask_bt + mask_life ;
     mask( isnan(mask) ) = 0 ;

% Generate Output variables

     % 1. column-indices of selected eddies (ie. eddy id's) : mask ~= 0

     eid        = find (mask) ; 
     estat.n    = size(eid,2) ;

     % 2. birth location of selected eddies

     ebx   = eddy.x(b_id) ; % first extract x-for birth position of all eddies
     ebx   = ebx(eid) ;     % extract again for selected eddies
 
     eby   = eddy.y(b_id) ;
     eby   = eby(eid) ;

     % 3. death location of selected eddies

     edx   = eddy.x(d_id) .* nod_byend ; % first extract x-for death position of all eddies
     edx   = edx(eid)                  ;     % extract again for selected eddies
 
     edy   = eddy.y(d_id) .* nod_byend ;
     edy   = edy(eid)                  ;

     % 4. life length/period of selected eddies

     elife = mask_life(eid) ;

     % 5. birth time of selected eddies

     ebt   = ref_li(eid) ; 

%-----
% DONE 
%-----

     return
