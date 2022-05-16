function  func_min_life ( grid, min_life, etr, tmpnc, outnc, bad)
%
%   Description : Apply minimum life cutoff to eddy tracking results.
%
%   Inputs : grid (struct)     - from func_init_grid
%            min_life (scalar) - (in No. of snapshots) minimum life cutoff for eddy
%            etr   (struct)    - information for writing global attributes
%            tmpnc (string)    - name of temporary eddy tracking NetCDF output file
%            outnc (string)    - name of final     eddy tracking NetCDF output file
%            bad   (scalar)    - missing_value field for NetCDF file
%            
%   Outputs : Trimmed outnc, created from tmpnc, after removing all eddies with
%                 lifetime < minimum_life.
%
%  Jaison Kurian
%  June/26/2008
%  
%  Aug/21/2009 : Added a new input grid, to facilitate the use of new 
%                      func_create_eddy_file which will write grid variables to
%                      output NetCDF file.  
%  Sep/24/2009 : Slightly more efficient identification of valid eddy ids.
%                      Few std printouts...for information on progress.
%-----------------------------------------------------------------------------
   global vnames
        
%------------
%  Initialize
%------------

    matmax = 1E25 ; % Set all variables values > this to NaN

    if ( nargin < 6 ) 
       error ([mfilename ':argchk'], '\n\t FATAL (%s) : Need 6 arguments, including grid struct.\n', mfilename)
    end  

%-----------------------------------------------------
%  find indices of eddies with life period >= min_life
%-----------------------------------------------------


     clear anti cyc

     fprintf ('\t NOTE (%s) : Reading ref var ....\n',mfilename) 

     ncp = netcdf ( tmpnc, 'read' ) ;

        anti_ref  = ncp{'anti_i'}(:) ;
        cyc_ref   = ncp{'cyc_i'}(:) ;

     ncp = close ( ncp ) ;

     fprintf ('\t NOTE (%s) : Finding Valid eddy IDs ....\n',mfilename) 

     anti_ref(abs(anti_ref)>matmax)    = 0 ;
     anti_ref(isnan(anti_ref))         = 0 ;
     cyc_ref(abs(cyc_ref)>matmax)      = 0 ;
     cyc_ref(isnan(cyc_ref))           = 0 ;

     cyc_long  =  find(sum(cyc_ref>0,2) >= min_life) ;
     anti_long =  find(sum(anti_ref>0,2) >= min_life) ;

     isz     = max ( length(anti_long), length(cyc_long) )  ;

     [nx nt]   = size(anti_ref) ;

%-----------------------------------------------------------------
%  Write eddies with life period >= min_life to final output file
%-----------------------------------------------------------------


     fprintf ('\t NOTE (%s) : Creating Destination NetCDF file ....\n',mfilename) 
     tmp1 = func_create_eddy_file (outnc, vnames, isz, etr, 0) ; 


     nvar   = length(vnames) ;

     ncp = netcdf ( tmpnc, 'read' ) ;
     nco = netcdf ( outnc, 'write' ) ;

        nco{'time'}(1:nt)      = ncp{'time'}(1:nt) ;
        nco{'xeddynum'}(1:isz) = [1:isz] ;
        
        tmp2 = sync ( nco ) ;

        for iv = 1:nvar ;
        fprintf ('\t NOTE (%s) : Writing  Var  %2i/%3i  %s ....\n',mfilename,iv, nvar, vnames{iv}) 
   
            tmp     =  ones(nt,isz) .* bad ;
            c_net   =  ['cyc_', vnames{iv}] ;
            tmpvar  = ncp{c_net}(:,:) ;
            tmpvar(tmpvar > matmax) = bad ;
            tmp(:,1:length(cyc_long)) = tmpvar(cyc_long,:)' ;
            nco{c_net}(1:nt,:) = tmp ;

            tmp     =  ones(nt,isz) .* bad ;
            a_net   =  ['anti_', vnames{iv}] ;
            tmpvar  = ncp{a_net}(:,:) ;
            tmpvar(tmpvar > matmax) = bad ;
            tmp(:,1:length(anti_long)) = tmpvar(anti_long,:)' ;
            nco{a_net}(1:nt,:) = tmp ;

        end  

     nco = close ( nco );
     ncp = close ( ncp ); 

    % delete (tmpnc) ;

     return
