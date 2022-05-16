function grid = func_init_grid_SanjivVersion (gridfile, reffile, omod)
%%
%% Description : Initialize ROMS RHO_grid, based on specified grid file.
%%                  No options for W-grid as of now.
%%
%% Input  : 1. grid file name, with full path (if not exist in pwd)
%%          2. reference ROMS output file, to retrieve theta_s, theta_b etc.
%%
%% Output : A structure, with the following variables:
%%
%%            grid.xkm     grid.sina
%%            grid.ykm     grid.cosa
%%            grid.area
%%            grid.lon
%%            grid.lat
%%            grid.depth
%%
%% Calls  : func_zlev
%%
%% UCLA ROMS Team
%% May/29/2008
%%
%%    Feb/04/2008 : added the option to get "basic" variables
%% Modification : Oct/04/2010 : f,pn, pm and areai are now part of basic package
%%-----------------------------------------------------------------------------

%---------------
%  Initial Check
%---------------

     omode = 'all' ;

     if ( nargin < 2 )
        error([mfilename ':argchk'], '\n   FATAL (%s) : Too few input arguments. \n', mfilename)
     elseif ( nargin == 3 )
        omode = lower(omod) ;
     end  

     if ( exist(gridfile,'file')  ~= 2 )
        error([mfilename ':argchk'], '\n   FATAL (%s) : Cannot find specified ROMS Grid file.\n \t\t %s \n', mfilename, gridfile)
     end

     if ( exist(reffile,'file')  ~= 2 )
        error([mfilename ':argchk'], '\n   FATAL (%s) : Cannot find specified ROMS output reference file.\n \t\t %s \n', mfilename, reffile)
     end

%-------------------------------------------------------------
%  Derive (rho_grid) : xkm, ykm, area, depth, cosa, sina, mask
%-------------------------------------------------------------

     m2km = 1/1000 ; % meter to km converter

     % do the basic part first
     ncgr = netcdf ( gridfile, 'read' ) ;
         grid.lon = ncgr{'lon_rho'}(:) ;
         grid.lat = ncgr{'lat_rho'}(:) ;
         grid.xkm = ncgr{'x_rho'}(:).* m2km ;
         grid.ykm = ncgr{'y_rho'}(:).* m2km ;
         grid.mask= ncgr{'mask_rho'}(:) ;
         grid.mask(grid.mask==0)  = NaN ;
         grid.pm       = ncgr{'pm'}(:) ; 
         grid.pn       = ncgr{'pn'}(:) ; 
         grid.f   = ncgr{'f'}(:) ;
         grid.areai      = grid.pm .* grid.pn ;

         if ( isempty(grid.xkm) || isempty(grid.ykm) || isempty(grid.mask) || isempty(grid.lon) || isempty(grid.lat) || isempty(grid.pn) || isempty(grid.pn) || isempty(grid.f) )
            error([mfilename ':readnc'], '\n   FATAL (%s) : Failed to retrieve variables from \n \t %s \n',mfilename, gridfile)
         end          

         if ( strcmp (omode, 'basic') ~= 1 )
       
             h        = ncgr{'h'}(:) ; 
             angle    = ncgr{'angle'}(:);
             mask     = ncgr{'mask_rho'}(:) ;
             masku    = ncgr{'mask_u'}(:) ;
             maskv    = ncgr{'mask_v'}(:) ;

             if ( isempty(grid.pm) || isempty(grid.pn) || isempty(h) || isempty(angle) || isempty(masku) || isempty(maskv) ) 
               error([mfilename ':readnc'], '\n   FATAL (%s) : Failed to retrieve variables from \n \t %s \n',mfilename, gridfile)
             end

             grid.cosa       = cos(angle); 
             grid.sina       = sin(angle);
             masku(masku==0) = NaN ;
             grid.masku      = masku ;      
             maskv(maskv==0) = NaN ;
             grid.maskv      = maskv ;      


             % coefficients to calculate grid depth
             ncrf = netcdf ( reffile, 'read' ) ;
                 Nsigma  = length(ncrf('s_rho')) ;
                 grid.N  = Nsigma ;
                 theta_s = ncrf.theta_s(:) ;
                 theta_b = ncrf.theta_b(:) ;
                 Tcline  = ncrf.hc(:)      ;
                 if ( isempty(Nsigma) || isempty(theta_s) || isempty(theta_b) || isempty(Tcline) )        
                    error([mfilename ':readnc'], '\n   FATAL (%s) : Failed to retrieve global attributes from \n \t %s \n',mfilename, reffile)
                 end
             ncrf = close (ncrf) ; 
             
             % compute depth 
             grid.depth = func_zlev(h,theta_s,theta_b,Tcline,Nsigma,'r');

          end
     ncgr = close (ncgr) ;  
     % clean up
     clear ncgr ncrf angle mask masku maskv m2km h  pm pn Nsigma theta_s theta_b Tcline ;
