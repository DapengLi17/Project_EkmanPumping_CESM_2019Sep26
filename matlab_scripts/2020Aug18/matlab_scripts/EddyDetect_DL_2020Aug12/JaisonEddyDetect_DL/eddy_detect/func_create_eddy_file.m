function [Missvalue vnames] = func_create_eddy_file ( fileout, vnames, dimsize, attr, flipdim )
%
%   function [Missvalue vnames] = func_create_eddy_file ( fileout, vnames, dimsize, attr, flipdim )
%
%   Description : Create eddy tracking output NetCDF file and write  
%                     header & coordinate variables.
%
%   Input : 1 fileout - (string, scalar)  output eddy tracking netcdf file
%           2 vnames  - (string, vector)  names of pre-defined output variables 
%                                         required
%           3 dimsize - (float,  scalar)  size of limited dimension (X or T)
%           4 attr    - (mixed, struct)   input parameters for eddy tracking.
%           5 flipdim - (float, scalar)   1 --> var(eid,time) with eid=unlimited 
%                                       0 --> var(time,eid) with time=unlimited 
%
%   Output : new fileout NetCDF file.
%            1. Missvalue --> missing value for NetCDF variables
%            2. vnames    --> updated vnames, if any of the basic ones are 
%                                     missing.
%
%   Modifications : Nov/29/2008 : A long awaited modification....option to choose 
%                     output variables through function argument "vnames".
%                     Few default variables will be added to make sure that the
%                     eddy tracking machienery will work fine.
%
%                     To add more variables, update vlist, ttl and unt, and
%                     vdef if required. Also update in func_get_eddies.m
%
%                   July/10/2009 : More generalizations to suit with new 
%                     eddy properties and tracking methods based on Q and SLA. 
%
%                     attr needs the following fields as struct :    
%                        tunits, torigin, tcalendar, method
%                   Oct/04/2010 : added new variables vort, vortmax, swirlmax
%____________________________________________________________________
% 
%----------------------------------------------------------------
   
     % check for input arguments

     if ( nargin < 5 )
          error('\n   FATAL (%s) : Too few input arguments (need 5). \n', mfilename)
     end     
 
     if ( isfield(attr,'method') == 1 )
          if ( strcmp(lower(attr.method),'q') == 1 )
                  unttmp = '1/s' ;
                %  if ( sum(strcmp(lower(vnames),'slaave')) == 1 || sum(strcmp(lower(vnames),'slamax')) == 1 )
                 %    error ('\n\t FATAL (%s) : Do not request SLA variables with Q-method.\n\n',mfilename)
                 % end
          elseif (strcmp(lower(attr.method),'sla') == 1 )
                  unttmp = 'cm' ;
          else 
              error ('\n\t  FATAL (%s) : Unknown eddy tracking method : %s. Only Q and SLA are valid methods.\n\n',mfilename, attr.method)
          end
     else 
         error ('\n\t  FATAL (%s) : Need to define eddy tracking method in attr struct.\n\n',mfilename)
     end 

     fprintf ( '\n \t  NOTE (%s) : Creating %s ......\n', mfilename, fileout)

     % Hard-wired time-axis part
   
     %tunits    = 'days since 1991-12-31 00:00:00' ; % moved to the area of 
     %torigin   = '31-DEC-1991 00:00:0' ;            % main script 
     %tcalendar = 'gregorian' ;
     Missvalue = -1.e+34 ; 

     % complete list of output variables, to check against specified ones.
     %   Their title and units to go with. All fields follow the same order.

     %vlist = {'i' 'j' 'x' 'y' 'xkm' 'ykm' 'r' 'rq' 'rkm' 'rkmq'  'd' 'dkm' 'dcokm', 'temp' 'salt' 'spice' 'vort' 'q' 'zeta' } ;
     vlist = {'i' 'j' 'x' 'y' 'xkm' 'ykm' 'r' 'rq' 'rkm' 'rkmq'  'd' 'dkm' 'dcokm'  'area' 'qave' 'qmax' 'slaave' 'slamax' 'cont' 'perr' 'vort' 'vortmax' 'swirlmax'} ;
%     vlist = {'i' 'j' 'x' 'y' 'xkm' 'ykm' 'r' 'rq' 'rkm' 'rkmq'  'd' 'dkm' 'dcokm'  'area' 'qave' 'qmax' 'slaave' 'slamax' 'sstave' 'sssave' 'cont' 'perr' 'vort' 'vortmax' 'swirlmax'} ;

     vdef   = {'i' 'j' 'x' 'y' 'xkm' 'ykm' 'r' 'rkm'} ; % default variables
     veav   = {'temp' 'salt' 'spice'  'q' 'zeta'}; % variable for eddy interioir averaging, not handled now
                                                         %   due to output filesize restrictions for long period tracking

     ttl{1}    = 'I Index of Nearest ROMS Grid point'                 ; % i
     ttl{2}    = 'J Index of Nearest ROMS Grid point'                 ; % j
     ttl{3}    = 'Center X position'                                  ; % x
     ttl{4}    = 'Center Y position'                                  ; % y
     ttl{5}    = 'Center X position in km'                            ; % xkm
     ttl{6}    = 'Center Y position in km'                            ; % ykm
     ttl{7}    = 'Radius (Vorticity Based)'                           ; % r
     ttl{8}    = 'Radius (Q Based)'                                   ; % rq        
     ttl{9}    = 'Radius (Vorticity Based) in km'                     ; % rkm
     ttl{10}   = 'Radius (Q Based) in km'                             ; % rkmq
     ttl{11}   = 'Distance traveled with respect to previous position'; % d
     ttl{12}   = 'Distance traveled with respect to previous position in km'; % dkm
     ttl{13}   = 'Distance from coast, along nearest J-line in km'    ; % dcokm
     ttl{14}   = 'Area covered by closed contour'                     ; % area
     ttl{15}   = 'Q-average within closed contour'                    ; % qave
     ttl{16}   = 'Q-Maximum within closed contour'                    ; % qmax
     ttl{17}   = 'SLA-average within closed contour'                  ; % slaave
     ttl{18}   = 'SLA-Maximum within closed contour'                  ; % slamax
     ttl{19}   = 'Contour value'                                      ; % slamax
     ttl{20}   = 'Percentage error in circular shape'                 ; % perr
   % ttl{14}   = 'Temperature (Eddy Interior Av.)'                    ; % temp   
   % ttl{15}   = 'Salinity (Eddy Interior Av.)'                       ; % salt 
   % ttl{16}   = 'Spiciness (Eddy Interior Av.)'                      ; % spice
   % ttl{18}   = 'Q-Parameter (Eddy Interior Av.)'                    ; % q   
   % ttl{19}   = 'Zeta (SSH) (Eddy Interior Av.)'                     ; % zeta
     ttl{21}   = 'Vorticity/f (Eddy Interior Av.)'                    ; % vort
     ttl{22}   = 'Max of Vorticity/f (Eddy Interior)'                 ; % vortmax
     ttl{23}   = 'Max of Swirl Velocity'                              ; % swirlmax

     unt{1}    = 'index'        ; % i
     unt{2}    = 'index'        ; % j
     unt{3}    = 'degrees_east' ; % x
     unt{4}    = 'degrees_north'; % y
     unt{5}    = 'km'           ; % xkm
     unt{6}    = 'km'           ; % ykm
     unt{7}    = 'Degrees'      ; % r
     unt{8}    = 'Degrees'      ; % rq   
     unt{9}    = 'km'           ; % rkm
     unt{10}   = 'km'           ; % rkmq
     unt{11}   = 'Degrees'      ; % d
     unt{12}   = 'km'           ; % dkm
     unt{13}   = 'km'           ; % dcokm
     unt{14}   = 'km^2'         ; % area
     unt{15}   = '1/s'          ; % qave
     unt{16}   = '1/s'          ; % qmax
     unt{17}   = 'cm'           ; % slaave
     unt{18}   = 'cm'           ; % slamax
     unt{19}   = unttmp         ; % cont
     unt{20}   = '%'            ; % perr
     unt{21}   = '1/s'          ; % vort
     unt{22}   = '1/s'          ; % vortmax
     unt{23}   = 'm/s'          ; % swirlmax
   % unt{14}   = 'Deg_C'        ; % temp 
   % unt{15}   = 'psu'          ; % salt 
   % unt{16}   = 'Sigma_Theta'  ; % spice
   % unt{17}   = '1/s'          ; % vort
   % unt{18}   = '1/s'          ; % q   
   % unt{19}   = 'm'            ; % zeta

     nlist     = length (vlist)  ; % number of variables in the list
     ndef      = length (vdef)   ; % number of default variables

        % always add default variables

     nvar      = length (vnames) ; % number of requested output variables
     for iv = 1:ndef 
        vex = [] ; 
        vex  = strmatch ( vdef{iv}, vnames, 'exact') ; % index of requested variable in vdef
        if ( isempty(vex) )
           nvar = nvar + 1 ;   
           vnames{nvar} = vdef{iv} ;
        end 
     end 

     % write header part of output file 

     nc = netcdf(fileout, 'clobber') ;
     if isempty(nc), error(' ## Could not open fileout for writing.'), end

         % dims
         xax  = 'xeddynum' ;
         time = 'time' ;
         
         % decide whether this call is for a tmp file or not
               % tmp file has dimension order (xeddynum, time)

         if (  flipdim == 1 )
	      nc(xax)  = 0 ;
              nc(time) = dimsize ;
              dim1     = xax ;
              dim2     = time ; 
         else 
	      nc(xax)  = dimsize ;
              nc(time) = 0 ;
              dim1     = time ;
              dim2     = xax ; 
         end
 
         % vars and atts
         nc{time}             = ncdouble(time) ;
         nc{time}.long_name   = ncchar('Time') ;
         nc{time}.units       = ncchar(attr.tunits) ;         
         nc{time}.time_origin = ncchar(attr.torigin) ;         
         nc{time}.calendar    = ncchar(attr.tcalendar) ;         
         nc{time}.axis        = ncchar('T') ;

         nc{xax}              = ncdouble(xax) ;
         nc{xax}.long_name    = ncchar('Eddy Identification Number') ;
	 nc{xax}.point_spacing= ncchar('even') ;
         nc{xax}.axis         = ncchar('X') ;

         vcount = 0 ;   % number of handled output variables
         for iv = 1:nvar 
            
            ivl = [] ;
            ivl = strmatch ( vnames{iv}, vlist, 'exact') ; % index of requested variable in list
            if ( ~isempty(ivl) && ivl <= nlist )

                c_net = ['cyc_', vnames{iv}] ;
                nc{c_net}               = ncfloat(dim1, dim2); 
                nc{c_net}.long_name     = ncchar ( [ 'CYC : ',  ttl{ivl} ] ) ;
                nc{c_net}.units         = ncchar ( unt{ivl} )                ;
                nc{c_net}.missing_value = ncfloat(Missvalue)                ;

                a_net = ['anti_', vnames{iv}] ;
                nc{a_net}               = ncfloat(dim1, dim2); 
                nc{a_net}.long_name     = ncchar ( [ 'ANTI : ', ttl{ivl} ] ) ;
                nc{a_net}.units         = ncchar ( unt{ivl} )                ;
                nc{a_net}.missing_value = ncfloat(Missvalue)                ;

                vcount = vcount + 1 ; 
            elseif ( ~isempty (strmatch ( vnames{iv}, veav, 'exact') ) ) ; % index of possible eddy interior averages
                fprintf ('\n\t WARNING (%s) : Eddy interior average should be \n\t\t found offline (var %s).\n', mfilename, vnames{iv})
   
            else
                fprintf ('\n\t WARNING (%s) : Requested variable %s is not valid and \n\t\t hence ignored.\n', mfilename, vnames{iv})

            end
         end
    
         % global atts

         wdir  = pwd;
         cdate = datestr(now, 'mmm/dd/yyyy HH:MM:SS AM') ;
         
         !hostname > hnametmp
         fidtmp   = fopen('hnametmp');
            hostname = fgetl(fidtmp);
         fclose(fidtmp);
         delete('hnametmp')
         
         nc.title       = ncchar(attr.title);
         nc.tools_used  = ncchar('Eddy Tracking');
         nc.computer = ncchar(hostname);
         nc.pwd      = ncchar(wdir);
         nc.grid     = ncchar(attr.grid) ;
         nc.datadir  = ncchar(attr.datadir) ;
         nc.inputfiles = ncchar([sprintf('%0.4d', attr.fids), ' - ', sprintf('%0.4d', attr.fide)]) ;
         nc.ztype    = ncchar(attr.ztype) ;
         nc.zvalue   = ncfloat (attr.zval) ;
         nc.contour_value_1E11 = ncfloat(attr.cont) ;
         nc.min_radius_km = ncfloat(attr.minrad);
         nc.mperr_area_percentage = ncfloat(attr.err);
         nc.nhann      = ncfloat(attr.nhann);
         if ( attr.nhann > 0 )
            smooth =   ['2D hanning, ', num2str(attr.nhann), ' times on vorticity, Q/W/SLA']  ;
         else 
            smooth = 'None';
         end
         nc.spatial_smoothing  = ncchar(smooth);
         nc.tracker_mfile = ncchar(attr.tracker) ;
         nc.method    =  ncchar(attr.method) ;
         nc.date     = ncchar(cdate);

         if ( flipdim == 1 ) 
              nc.message  = ncchar('#######   This is a TEMPORARY eddy tracking file   ##########');
         end

     nc = endef(nc) ; % leave define mode

     nc = close (nc) ; % close netcdf file

%-----
% DONE            
%-----
