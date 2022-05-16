function grid = func_init_grid (gridfile, reffile, omod)
%
% function grid = func_init_grid (gridfile, reffile, omod)
%
% Description : Initialize ROMS RHO_grid, based on specified grid file.
%                  No options for W-grid as of now.
%
% Input  : 1. grid file name, with full path (if not exist in pwd)
%          2. reference ROMS output file, to retrieve theta_s, theta_b etc.
%          3. grid initialization mode ('basic' or 'all')
%
% Output : A structure, with the following variables:
%
%            grid.xkm     grid.sina
%            grid.ykm     grid.cosa
%            grid.area
%            grid.lon
%            grid.lat
%            grid.depth
%
% Calls  : func_zlev

%_______________________________________________________________________
%  Written By     : Jaison Kurian (jaisonk@tamu.edu
%                                  jaisonkuriann@gmail.com) 
%                   Oceanography Dept., TAMU
%  Date           : Jan/11/2011
%  eTRACK Version : 2.00 
%  Ref. Version   : etrack_07_sla_circ
%  Copyright      : See etrack_manual.pdf
%  Modifications  : see below
%_______________________________________________________________________
%  Feb/04/2008 : added the option to get "basic" variables
%  Oct/04/2010 : f,pn, pm and areai are now part of basic package
%  Jul/09/2014 : Following fields are computed if missing: pn,pm,x_rho,y_rho,f.
%_______________________________________________________________________

%---------------
%  Initialize
%---------------

     %-----------nargin chk---------------------------------------------%
     omode = 'all' ;
     if ( nargin < 2 )
        error('\n   FATAL (%s) : Too few input arguments. \n', mfilename)
     elseif ( nargin == 3 )
        omode = lower(omod) ;
     end  

     %-----------file existence chk-------------------------------------%
     if ( exist(gridfile,'file')  ~= 2 )
        error('\n   FATAL (%s) : Cannot find specified ROMS Grid file.\n \t\t %s \n', mfilename, gridfile)
     end
     if ( exist(reffile,'file')  ~= 2 )
        error('\n   FATAL (%s) : Cannot find specified ROMS output reference file.\n \t\t %s \n', mfilename, reffile)
     end

%-------------------------------------------------------------
%  Derive (rho_grid) : xkm, ykm, area, depth, cosa, sina, mask
%-------------------------------------------------------------

     m2km = 1/1000 ; % meter to km converter

     %-----------Open and Read Grid file--------------------------------%
     ncgr = netcdf ( gridfile, 'read' ) ;
         grid.lon      = ncgr{'lon_rho'}(:) ;
         grid.lat      = ncgr{'lat_rho'}(:) ;
         grid.xkm      = ncgr{'x_rho'}(:).* m2km ;
         grid.ykm      = ncgr{'y_rho'}(:).* m2km ;
         grid.mask     = ncgr{'mask_rho'}(:) ;
         grid.mask(grid.mask==0)  = NaN ;
         grid.pm       = ncgr{'pm'}(:) ; 
         grid.pn       = ncgr{'pn'}(:) ; 
         grid.f        = ncgr{'f'}(:) ;
         grid.areai    = grid.pm .* grid.pn ;

         if ( isempty(grid.mask) || isempty(grid.lon) || isempty(grid.lat) || isempty(grid.pn) || isempty(grid.pm) || isempty(grid.f) )
%            error('\n   FATAL (%s) : Failed to retrieve variables from \n \t %s \n',mfilename, gridfile)
         end          
     %-----------Find pm and pn if needed-------------------------------%
                      % step 1. use great-circle-distance method from func_gcdist.m
                      % step 2. convert km to meters
                      % step 3. values at at edges equal to that at the closest interior point
         if ( isempty(grid.pm) || isempty(grid.pn) )

             fprintf('\n   WARNING (%s) : Missing pm & pn being estimated from lon & lat.\n',mfilename)

             km2m        = 1000.0     ; % km to meter converter
             [jsz isz]   = size(grid.lon) ;

             dx(:,2:isz) = func_gcdist(grid.lat(:,2:isz),grid.lon(:,2:isz),grid.lat(:,1:isz-1),grid.lon(:,1:isz-1)) ; % step 1
             dx          = dx .* km2m   ; % step 2 
             dx(:,1)     = dx(:,2)      ; % step 3 
             dx(:,isz)   = dx(:,isz-1)  ; % step 3
             grid.pm     = 1./dx        ;

             dy(2:jsz,:) = func_gcdist(grid.lat(2:jsz,:),grid.lon(2:jsz,:),grid.lat(1:jsz-1,:),grid.lon(1:jsz-1,:)) ;
             dy          = dy .* km2m   ; % step 2
             dy(1,:)     = dy(2,:)      ; % step 3
             dy(jsz,:)   = dy(jsz-1,:)  ; % step 3
             grid.pn     = 1./dy        ;

         end
     %-----------Find x_rho and y_rho (xkm & ykm) if needed-------------%
         if ( isempty(grid.xkm) || isempty(grid.ykm) )
             fprintf('\n   WARNING (%s) : Missing x_rho & y_rho being estimated from pm & pn.\n',mfilename)
             grid.xkm          = cumsum(1./grid.pm,2) .* m2km ;
             grid.xkm(:,2:end) = grid.xkm(:,1:end-1) ;
             grid.xkm(:,1)     = 0 ;  
            
             grid.ykm          = cumsum(1./grid.pn,1) .* m2km ;
             grid.ykm(2:end,:) = grid.ykm(1:end-1,:) ;
             grid.ykm(1,:)     = 0 ;
         end
      
     %-----------Find f if needed---------------------------------------%
                       % Step 1. f = 2 * sin (lat .* deg2rad);  OMEGA=7.29E-5
                       % Ref : http://en.wikipedia.org/wiki/Coriolis_effect#Formula 
         if ( isempty(grid.f) )
             fprintf('\n   WARNING (%s) : Missing f being estimated from lat & lon.\n',mfilename)
             deg2rad    = pi/180             ; % radians
             grid.f     = 2 .* 7.29E-5 .* sin(grid.lat .* deg2rad); 
         end

     %-----------Mode check---------------------------------------------%
         if ( strcmp (omode, 'basic') ~= 1 )
       
             h        = ncgr{'h'}(:) ; 
             angle    = ncgr{'angle'}(:);
             mask     = ncgr{'mask_rho'}(:) ;
             masku    = ncgr{'mask_u'}(:) ;
             maskv    = ncgr{'mask_v'}(:) ;

             if ( isempty(grid.pm) || isempty(grid.pn) || isempty(h) || isempty(angle) || isempty(masku) || isempty(maskv) ) 
               error('\n   FATAL (%s) : Failed to retrieve variables from \n \t %s \n',mfilename, gridfile)
             end

             grid.cosa       = cos(angle); 
             grid.sina       = sin(angle);
             masku(masku==0) = NaN ;
             grid.masku      = masku ;      
             maskv(maskv==0) = NaN ;
             grid.maskv      = maskv ;      

     %-----------Open/Read/Close reffile--------------------------------%
             ncrf = netcdf ( reffile, 'read' ) 
                 Nsigma  = length(ncrf('s_rho')); 
                 grid.N  = Nsigma ;
                 theta_s = ncrf{'theta_s'}(:) ;
                 theta_b = ncrf{'theta_b'}(:) ;
                 Tcline  = ncrf{'hc'}(:)      ;

                 if ( isempty(Nsigma) || isempty(theta_s) || isempty(theta_b) || isempty(Tcline) )        
                    error('\n   FATAL (%s) : Failed to retrieve global attributes from \n \t %s \n',mfilename, reffile)
                 end
             ncrf = close (ncrf) ; 
             
     %-----------Compute Depth------------------------------------------%
             grid.depth = func_zlev(h,theta_s,theta_b,Tcline,Nsigma,'r');

          end
     %-----------Close grid file----------------------------------------%
     ncgr = close (ncgr) ;  

     %-----------Cleanup------------------------------------------------%
     clear ncgr ncrf angle mask masku maskv m2km h pm pn Nsigma theta_s theta_b Tcline dx dy ;

%--------
% DONE
%--------
     return
