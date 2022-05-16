function diag = func_get_diag_sla (grid, prog, climfile, vnames, nhann)
%
% Description : Matlab function to calculate model diagnostics required for eddy
%                 tracking.
%
% Input :  1. grid - struct, ROMS grid fields, result from calling func_init_grid.
%          2. prog - struct, ROMS prognostic fields, result from calling 
%                         func_get_hslice (variables on their native grid)
%          3. climfile - string, Seasonal climatology file
%          4. vnames   - string, 1D-vector, names of required output variable
%          5. nhann    - scalar, float, apply hanning 2D smoother this many
%                            times
%
% Output : diag, struct, with following fields (matrix row and column length
%                will be 1 less than that of input)
%             diag.vort : vorticity (1/s), using backward finite differencing.
%             diag.q    : Q-parameter described by Isren-Fontanet et al (JAOT, 2003).
%             diag.w    : W-paremeter described by Chelton et al., (GRL, 2007)
%             diag.sla  : SLA, with seasonal mean removed
%
% Ref. Isren-Fonanet et al., 2003, Identification of marine eddies from altimetric
%               maps, JAOT, Vol.20, pp.772-778.
%      Chelton et al., 2007, Global observations of large oceanic eddies, GRL,
%               Vol.34,L15606,doi:10.1029/2007GL030812
%      Jackett and McDougall, 1995, Minimal Adjustment of Hydrographic Profiles to 
%               Achieve Static Stability, Vol.12(4), pp.381???389.            
%      Flament, 2002, A state variable for characterizing water masses and their 
%               diffusive stability: spiciness, Prog. Oceanogr. Vol.54(1-4,
%               July-Sep), pp.493-501.
%
% Written By : ROMS Team, UCLA
% Written On : May/31/2008
%
% Modifications : July/30/2008, hanning smoothing is standard now. Added a new
%                               input argument "nhann".
%                 Jan/20/2009,  added options to output normalized SLA fields.
%                 Feb/04/2009, basic version, ROMS SLA, no normalization or av.
%                                removal
%                 April/14/2009 added few standard diagnostic variables for 
%                                 purposes with the Q-based method.
%
%                 July/01/2009  added vnames
%                 Nov/10/2014   Added option for "no climatology file".
%-------------------------------------------------------------------------------

%-----------
% Initialize 
%-----------


     if ( nargin < 4)
        error([mfilename ':argchk'], '\n   FATAL (%s) : Too few input arguments.\n', mfilename)
     end

     if ( nhann < 0 ) 
        error([mfilename ':argchk'], '\n   FATAL (%s) : nhann should be >= 0.\n', mfilename)
     end

     if ( ~isempty(climfile) )
        if ( exist(climfile,'file')  ~= 2 )
            error([mfilename ':argchk'], '\n   FATAL (%s) : Cannot find specified ROMS Clim file.\n \t\t %s \n', mfilename, climfile)
        end
     end

     vnames = lower (vnames) ;

%-----------------
% Find diagnostics
%-----------------
     % SSH
     diag.sla = prog.zeta * 100 ; % m to cm

    % SST
%      diag.sst = prog.sst ;

    % SSS
%      diag.sss = prog.sss; 

    % UTAU, VTAU, QNET
%      diag.utau= prog.utau;
%      diag.vtau= prog.vtau;
%      diag.qnet= prog.qnet;

    % SBL
%      diag.sbl = prog.sbl; 

    % Mean SSH

     if ( isempty(climfile) )
        ssh = 0.0 ;
     else 
        nc = netcdf ( climfile, 'read' ) ;
           ssh = squeeze ( nc{'zeta'}(prog.ssn,:,:) ) .* grid.mask * 100; % m to cm
        nc = close ( nc ) ;   
     end
     
     % SLA

     diag.sla  = (diag.sla - ssh) ; % remove time mean sea surface height (SSH)
     
     need_der = false ;

     vdiag = {'q' 'w' 'vort'} ;
     for iv = 1:length(vnames) 
        if ( sum(strcmp(vnames,vdiag{iv})) == 1 ) ; need_der = true ; end
     end

     if ( need_der )
        % vorticity, q & w
        
        [jsz isz] = size(prog.zeta) ;

        [geo.u geo.v] = func_roms_geovel (grid.f, diag.sla./100, grid.pm, grid.pn ) ;
        
        udx  = 2 * geo.u ./ ( grid.pm(:,1:isz-1) + grid.pm(:,2:isz) ) ; % u * av dx  -> m^2/s
        vdx  = 2 * geo.v ./ ( grid.pm(1:jsz-1,:) + grid.pm(2:jsz,:) ) ; % v * av dx  -> m^2/s       vom
        udy  = 2 * geo.u ./ ( grid.pn(:,1:isz-1) + grid.pn(:,2:isz) ) ; % u * av dy  ->             uon
        vdy  = 2 * geo.v ./ ( grid.pn(1:jsz-1,:) + grid.pn(2:jsz,:) ) ; % v * av dy  -> m^2/s

        
        ddx_u = zeros ( jsz, isz ) ; 
        ddx_u(2:end-1,2:end-1) = grid.areai(2:end-1,2:end-1) .* ( udy(2:end-1,2:end)-udy(2:end-1,1:end-1) ) ;
        ddx_v = grid.areai .* psi2rho ( (vdy(:,2:isz) - vdy(:,1:isz-1)) ) ;
        ddy_u = grid.areai .* psi2rho ( (udx(2:jsz,:) - udx(1:jsz-1,:)) ) ;

        if ( sum(strcmp(vnames,'vort')) == 1 )
           diag.vort  = ( ddx_v - ddy_u ) ./ grid.f ;          % vorticity normalized by f
        end
        if ( sum(strcmp(vnames,'q')) == 1 )
           diag.q     = 1E11 * ( -1 .* ddx_u.^2 - (ddx_v .* ddy_u) ) ;
        end 
        if ( sum(strcmp(vnames,'q')) == 1 )
           diag.w     = 1E11 * 4 * ( ddx_u.^2 + (ddx_v .* ddy_u) ) ;
        end

     end

     if ( sum(strcmp(vnames,'spd')) == 1 )
         diag.spd  = ( u2rho_2d(geo.u).^2 + v2rho_2d(geo.v).^2 ).^0.5 ;
     end

     if ( nhann > 0 ) % Smoothing
       for iv = 1:length(vnames)
          vnam   = ['diag.' vnames{iv}] ;  
          for nh = 1:nhann
             eval( [ vnam, '  = func_hann2d_fast(', vnam, ')', ';'] ) ;
          end
       end
     end
  

% DONE

     return 
