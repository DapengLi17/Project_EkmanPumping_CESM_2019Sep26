function eddy = func_init_eddy ( etype, efile, vnames, tids, tide )
%
%  Description : Matlab function to initialize eddy tracking results for a
%                   given eddy type.
%
%  Inputs  : 1. etype - (string, scalar) eddy type, 'cyc' or 'anti'
%            2. efile - (string, scalar) eddy tracking output NetCDF filename
%            3. vnames- (string, vector) variables to be initialized from eddy
%                                          tracking output file.
%            4. tids  - (float, scalar, optional)  starting time index to 
%                                extract data. If not specified, then set to 1.
%            5. tide  - (float, scalar, optional)  ending   time index to
%                                extract data. If not specified, then set to max.
%    
%  Outputs : 1. eddy  - (float, struct)  eddy variables requested via vnames,
%                             with proper NaNs and for the specified time period.
%
%                        Requested or not, "i" (cyc_i or anti_i) will be included
%                        in the output struct, to make life easier.
%

%_____________________________________________________________________________
% This is a part of Eddy Tracking Tool Kit, Copyright @ UCLA ROMS Team, 2008
% Written By : Jaison Kurian (jaison@atmos.ucla.edu)
% Written On : Nov/05/2008
%_____________________________________________________________________________

 
% Initialize

     if ( nargin < 3 || nargin == 4)
        error([mfilename ':argchk'], '\n   FATAL (%s) : Wrong number of input arguments. Need either 3 or 5.\n', mfilename)
     end  

     if ( ~ischar(etype) || ~ischar(efile) ) 
        error([mfilename ':argchk'], '\n   FATAL (%s) : Wrong input type (arg1/arg2) - Need strings, \n\t  eddy type and eddy filename.\n', mfilename)
     end

     if ( ~strcmp(etype,'cyc') && ~strcmp(etype,'anti') )
        error([mfilename ':argchk'], '\n   FATAL (%s) : Wrong eddy type (%s, arg3) \n\t\t - Need either cyc or anti.\n', mfilename, etype)
     end

     if ( ~ischar(vnames{1}) )
        error([mfilename ':argchk'], '\n   FATAL (%s) : Wrong input type (arg3) - need strings as a vector, \n\t list of eddy tracking variables to be initialized.\n', mfilename)
     end  

     if ( exist(efile,'file')  ~= 2 )
        error([mfilename ':argchk'], '\n   FATAL (%s) : Cannot find specified eddy tracking output file.\n \t\t %s \n', mfilename, efile)
     end   

     tsize = [] ;
     nce = netcdf ( efile, 'read' ) ;
         tsize = size ( nce{'anti_i'}(:), 1) ; 
     nce = close ( nce ) ;
     if ( isempty(tsize) )
        error([mfilename ':argchk'], '\n   FATAL (%s) : Cannot extract info on anti_i from\n \t\t %s \n', mfilename, efile)
     end

     if ( nargin > 3 )
        if ( ~isfloat(tids) || ~isfloat(tide) || ~func_shape(tids,0) || ~func_shape(tide,0) )
           error([mfilename ':argchk'], '\n   FATAL (%s) : Wrong arg type (args 4 & 5) : Need float scalars, \n\t\t starting and ending time indices).\n', mfilename)
        end 
        if ( min(tids,tide) < 1 || max(tids,tide) > tsize )
           error([mfilename ':argchk'], '\n   FATAL (%s) : Wrong time range (args 4 & 5) : should be \n\t\t between 1 and %g (present range %g:%g).\n', mfilename, tsize, tids, tide)
        end
     else 
        tids = 1 ;
        tide = tsize ;  
     end

% Initialize requested eddy variables

     matmax = 1.e+15 ;
     nvar = length(vnames) ;

   %  iex  = sum(strcmp (vnames, 'i')) ;
   %  if ( iex < 1 )
   %     nvar         = nvar + 1 ;
   %     vnames{nvar} = 'i' ;
   %  end
 
     % loop through requested variables     
     nce = netcdf ( efile, 'read' ) ;

        for iv = 1:nvar
           tmpvar = [] ;  bad = [] ;
           v_net  = [etype, '_', vnames{iv}] ;
           v_mat  = ['eddy.', vnames{iv}] ;
           tmpvar = squeeze (nce{v_net}(tids:tide,:)) ;
           if ( isempty (tmpvar) )
              error([mfilename ':ncvarchk'], '\n   FATAL (%s) : Cannot retrieve %s from \n\t\t %s.\n', mfilename, v_net, efile)
           end
           bad    = nce{v_net}.missing_value(:) ;
           if ( isempty (bad) )
              bad = -1.e+34;
           end
           tmpvar(tmpvar == bad) = NaN    ; % this step is not necessary due to the
           tmpvar(abs(tmpvar) > matmax) = NaN ; % the check with matmax
           eval( [v_mat, ' = tmpvar ;' ]) ;
        end
        
     nce = close ( nce ) ;
     
     % get eddy type to a field, for later use
     v_mat = ['eddy.', 'type'] ;  
     eval( [v_mat, ' = etype ;' ]) ;

% DONE
    
     return
