function  hslice = func_get_hslice (grid, filein, tindx, ztype, zval, dbm)
%
% Description : Function to cut an XY-slice, at given time & depth/sigma/
%                 isopycnal (rather sigma_theta), from 4D (XYZT) roms output.
%                 Output is on native grid rho/u/v grid itself.
%
% Input  : 1. grid   - struct, created by calling func_init_grid
%          2. filein - string, input ROMS avg filename
%          3. tindx  - scalar, time index with respect to filein
%          4. ztype  - string : 'depth' --> vertical mode is depth  OR
%                               'sigma' --> vertical mode is sigma layer number  OR
%                               'isopyc' --> vertical mode is isopycnal (sigma_theta
%                               'none' or 'zeta' or 'ssh' --> to extract just zeta.
%          5. zval   - scalar : meters (-ve) if ztype is 'depth'
%                               sigma layer number (+ve) if ztype is 'sigma' 
%                               isopycnal value (+ve) if ztype is 'isopyc' 
%          6. dbm    - string : to override global field 'dbmod' (see below)
%   
%
% Global : dbmod (scalar) : debug mode, 1 for debugging and 0 for not
%                                                              
% Output : struct hslice, with following fields (mask applied)
%            hslice.temp, hslice.salt, hslice.zeta
%            hslice.u & hslice.v   (both on rho grid)
%
% Needs  : func_zsurf.m  func_densityjm95 u2rho_3d v2rho_3d
%
% Written By : UCLA ROMS Team (jaison@atmos.ucla.edu)
% Written On : May/29/2008 
% Copyright @ UCLA ROMS Team
% Tool       : Eddy Tracker
% Version    : 1.0
%
% July/01/2009 : If ztype is 'none' or 'ssh' or 'zeta', then only the variable 
%                      'zeta' is returned.
%                              
%---------------------------------------------------------------------------------
    global dbmod
%-----------
% Initialize 
%-----------

     if ( ~exist('dbm', 'var') ) ; dbm = dbmod ; end  % choose debugging mode
     vaxis = lower ( ztype ) ;

     if ( (strcmp(vaxis,'none') == 1 || strcmp(vaxis,'zeta') == 1) || strcmp(vaxis,'ssh') == 1)
          vaxis = 'zeta' ;
          zval  = 0 ;
          zsize = grid.N ;
     else 
          zsize = grid.N ;
     end


     if ( dbm == 1 ) % if debugging mode on
        if ( nargin < 5 )
           error([mfilename ':argchk'], '\n   FATAL (%s) : Too few input arguments. \n', mfilename)
        end  

        if ( exist(filein,'file')  ~= 2 )
           error([mfilename ':argchk'], '\n   FATAL (%s) : Cannot find specified ROMS Avg file.\n \t\t %s \n', mfilename, filein)
        end     

        ncr = netcdf ( filein, 'read' ) ;
           tsize = [] ; 
           tmp1  = length( ncr{'ocean_time'}(:) ) ;  tmp2 = length( ncr{'scrum_time'}(:) ) ;
           tsize = max (tmp1, tmp2) ; 
           if ( tsize < 1 )
              error([mfilename ':timechk'], '\n   FATAL (%s) : Cannot get valid time axis size from specified ROMS Avg file.\n \t\t %s \n', mfilename, filein)
           end
        ncr = close (ncr) ;
     
        if ( tindx < 1 || tindx > tsize ) 
            error([mfilename ':timechk'], '\n   FATAL (%s) : Requested time index (%g) is outsize actual range (1 - %g) in \n \t\t %s \n', mfilename, tindx, tsize, filein)
        end 

        switch vaxis
           case {'depth' 'dep'}
               zmax = max(max(max( abs( grid.depth ) ))) ;
               if ( abs(zval) > zmax )
                  error([mfilename ':zchk'], '\n   FATAL (%s) : %g is unrealistic compared to Maximum depth of current grid (%g).\n', mfilename, zval, -1*zmax)
               end
           case {'sigma' 'sig'}
               if ( zval < 1 || zval > zsize || floor(zval) ~= zval )
                   error([mfilename ':zchk'], '\n   FATAL (%s) : %g is not a valid entry for Sigma Z-level (range is 1 - %g) \n', mfilename, zval, zsize)
                end
           case {'isopy' 'iso' 'isopyc' 'isopycnal'}
               if ( zval < 20 || zval > 30 )
                   error([mfilename ':zchk'], '\n   FATAL (%s) : %g is not a valid entry for Isopycnal Z-surface (range is 20 - 30) \n', mfilename, zval)
                end
           case {'zeta' 'ssh' 'none'} 
                % do nothing....just go to end
           otherwise 
               error([mfilename ':ztypechk'], '\n   FATAL (%s) : %s is not a valid zsurface type (use depth, sigma, or isopyc). \n', mfilename, zval)
        end

     end 

%--------------
% Preprocessing
%--------------

     ncr = netcdf ( filein, 'read' );
        if ( zval ~= 0 )  % 
           u_3d  = squeeze( ncr{'u'}(tindx,:,:,:) ) ;  
           v_3d  = squeeze( ncr{'v'}(tindx,:,:,:) ) ;  
           te_3d = squeeze( ncr{'temp'}(tindx,:,:,:) ) ;
           sa_3d = squeeze( ncr{'salt'}(tindx,:,:,:) ) ;
           rsize = size(te_3d) ; 
        end
        hslice.zeta = squeeze( ncr{'zeta'}(tindx,:,:) ) .* grid.mask  ; % it is a 2D var
        hslice.sst  = squeeze( ncr{'temp'}(tindx,zsize,:,:) ) .* grid.mask  ; % it is a 2D var
        hslice.sss  = squeeze( ncr{'salt'}(tindx,zsize,:,:) ) .* grid.mask  ; % it is a 2D var
        hslice.utau = u2rho_2d( squeeze( ncr{'sustr'}(tindx,:,:) ) ) .* grid.mask  ; % it is a 2D var
        hslice.vtau = v2rho_2d( squeeze( ncr{'svstr'}(tindx,:,:) ) ) .* grid.mask  ; % it is a 2D var
        hslice.qnet = squeeze( ncr{'shflux'}(tindx,:,:) ) .* grid.mask  ; % it is a 2D var

        hslice.sbl  = squeeze( ncr{'Hsbl'}(tindx,:,:) ) .* grid.mask  ; % it is a 2D var


     ncr = close (ncr) ;  

%------------------------------------------------
% Select Vertical axis type and extract variables
%------------------------------------------------


     switch vaxis

         case {'depth' 'dep'}

           zval = abs (zval) * -1  ;  % make it -ve (is it required ??)
           depth_u     = rho2u_3d ( grid.depth ) ;
           depth_v     = rho2v_3d ( grid.depth ) ;
           hslice.u    = func_zsurf ( u_3d, depth_u, zval ) .* grid.masku ;
           hslice.v    = func_zsurf ( v_3d, depth_v, zval ) .* grid.maskv ;
           hslice.temp = func_zsurf ( te_3d, grid.depth, zval ) .* grid.mask ;
           hslice.salt = func_zsurf ( sa_3d, grid.depth, zval ) .* grid.mask ;

        case {'sigma' 'sig'}  

           hslice.u    = squeeze ( u_3d(zval,:,:) ) .* grid.masku ;
           hslice.v    = squeeze ( v_3d(zval,:,:) ) .* grid.maskv;
           hslice.temp = squeeze ( te_3d(zval,:,:) ) .* grid.mask ;
           hslice.salt = squeeze ( sa_3d(zval,:,:) ) .* grid.mask ;

        case {'isopy' 'iso' 'isopyc' 'isopycnal'} 

           sig_3d      = func_densityjm95 ( sa_3d, te_3d, 0 ) - 1000.0 ;
           sig_3du     = rho2u_3d  ( sig_3d ) ;
           sig_3dv     = rho2v_3d  ( sig_3d ) ;
           hslice.u    = func_zsurf ( u_3d, sig_3du, zval ) .* grid.masku ;
           hslice.v    = func_zsurf ( v_3d, sig_3dv, zval ) .* grid.maskv ;
           hslice.temp = func_zsurf ( te_3d, sig_3d, zval ) .* grid.mask ;
           hslice.salt = func_zsurf ( sa_3d, sig_3d, zval ) .* grid.mask ;

     end    
     

