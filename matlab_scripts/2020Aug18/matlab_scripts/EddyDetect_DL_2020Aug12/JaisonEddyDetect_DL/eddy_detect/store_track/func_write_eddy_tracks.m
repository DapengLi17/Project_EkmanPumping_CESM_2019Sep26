function func_write_eddy_tracks (cyctrack, antitrack, tracknc, tmpnc, vnames, bad)
%
%  function func_write_eddy_tracks (cyctrack, antitrack, tracknc, tmpnc, vnames, bad)
%
%  Description : Matlab function to write eddy tracking results to tmp NetCDF file.
%                    Uses new version of func_track_eddy (2009/01/02).
%  
%  Inputs : cyctrack  (float, struct) : output from func_track_eddy, for cyc
%           antitrack (float, struct) : output from func_track_eddy, for anti
%           tracknc   (character, scalar) : name of output tmp track NetCDF file
%           tmpnc     (character, scalar) : name of tmp eddy detection NetCDF file
%           vnames    (character, 1D vector) : list of eddy variables 
%           bad       (float, scalar)        : NetCDF missing_value flag
%
%  Output : Updated values for vnames in tracknc
%
%  Jan/05/2008
%  Jan/06/2008 : Need dkm in vnames. Added writing for dkm from tracking outputs.

    
    if ( sum(strcmp (vnames, 'dkm')) ~= 1 ) ;
         error([mfilename ':argchk'], '\n   FATAL (%s) : Cannot find variable dkm in vnames.\n', mfilename)
    else 
         vids = find(~strcmp(vnames,'dkm')) ; % if dkm is present in the vnames
         vnames = vnames(vids) ;              %   remove it, since it is added later  
    end


    etypes = {'cyc' 'anti'} ;
    nvars = length (vnames) ;

    nc   = netcdf ( tmpnc, 'read' ) ;
    nctr = netcdf ( tracknc, 'write' ) ;

       tsize   = cyctrack.tsz ;
       ntracks = max(cyctrack.n, antitrack.n);
      
       nctr{'xeddynum'}(1:ntracks) = [1:1:ntracks] ;
       nctr{'time'}(1:tsize)       = nc{'time'}(:) ;
       err = sync (nctr) ;
      
       for ie = 1:2  % iterate through eddy types
      
          etype   = etypes{ie} ;
          isize   = eval([etype 'track.isz']) ;
          ntr     = eval([etype 'track.n']) ;

          for itr = 1:ntr                                % iterate through valid tracks
              ts   = eval([etype 'track.ts(itr)'])  ;
              te   = eval([etype 'track.te(itr)'])  ;
              lind = eval([etype 'track.lind(itr,1:te-ts+1)'])  ;
              for iv = 1:nvars                           % iterate through each variable name
                tmp        = zeros(1,tsize) + bad      ; % init with bad values
                net_v      = [etype '_' vnames{iv}]    ; % extract data from tmp file 
                tmp_all    = nc{net_v}(1:isize,:)      ; % account for re-scaling

                tmp(ts:te) = tmp_all (lind)            ; 
                nctr{net_v}(itr,:) = tmp               ;
              end                                        % end iteration through variable names
              tmp        =  zeros(1,tsize) + bad      ; % init with bad values
              tmp(ts:te) =  eval([etype 'track.dkm(itr,1:te-ts+1)'])  ;
              nctr{[etype '_dkm']}(itr,:) = tmp        ;
          end                                            % end iteration through tracks
      
          if ( ntr < ntracks )                           % for "trailing dummy" tracks 
             tmp = zeros(1,tsize) + bad                ; % write bad values
             for itr = ntr+1:ntracks
                for iv = 1:nvars
                   net_v      = [etype '_' vnames{iv}] ; 
                   nctr{net_v}(itr,:) = tmp            ;
                end
                nctr{[etype '_dkm']}(itr,:) = tmp ; 
             end
          end
      
       end % end eddy type

% DONE 

    nctr = close (nctr) ; 
    nc   = close (nc) ;

   
