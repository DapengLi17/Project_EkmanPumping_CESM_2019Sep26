%
%  Description : Eddy tracker, SLA-based method for ROMS.
%
%                  Please see the README file for the details of this version.
%
%  Updated : July/01/2009
%              Oct/06/2010 : see eddy_viewer for updates.
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    clear all ; %close all  
%---------------------------------------------------------------------------------------
    global omode dbmod vnames 
%---------------------------------------------------------------------------------------

    omode       = 'fast' ; 
    dbmod       = 0 ;


    contour_vals= [-50:1:50] ; % contour values of normalized (by std of SSH) SLA (cm)
                                  %     to search for
    min_radius  = 45 ; %10      ; % in km , required minimum radius to define an eddy (>=)
    max_radius  = 150 ; %200      ; % in km , required maximum radius to define an eddy (>=)
    min_points  = 1 ;
    mperr_area  = 50       ; % max percentage error in area, for circle fitting to Q
    nhann       = 0        ; % apply 2D hanning smoother on vort., Q/W this many times
    min_slaamp  = 1        ; % minimum SLA amplitude in cm
    dtime       = 2 ; % spacing between snapshots 
    min_life    = 1        ; % in number of snapshots, minimum life period for eddy
    leave       = 2        ; % for joining brocken tracks

    zvalue      = 0      ; % z-axis/surface value, according to ztype
    ztype       = 'none' ; % depth/sigma/isopyc/none/zeta/ssh


    romsgrid    = '/home/jaison/work/store/usw51_grd.nc' ; % ROMS Grid file (for model domain)

    directory   = '/mnt/paracas/capet_taniwha/USW5.1/CLIM_QSCATXA_SODA_HFLUX' ;
    clim_file   = '/home/jaison/work/store/usw51_Smean.nc' ; % Sclim

    model       = 'usw51'   ; % prefix for input files, typically model name
    filetype    = 'avg'     ; % average/history
    fid_start   =  90       ; % id for first file
    fid_end     =  1605 ; %2970     ;


    tinit       = [1901 1 1] ; % model start time, for finding seasons
    tcalendar   = '360_DAY' ;                      % 
 
    vnames      = {'r' 'x' 'y' 'i' 'j' 'rkm' 'xkm' 'dkm' 'ykm' 'dcokm' 'area' 'slaave' 'slamax' 'cont' 'perr' 'vort' 'vortmax' 'swirlmax'} ;
      
    
    fileout     = 'etrack_usw51_sla.nc' ;
    filemerge   = 'etrack_usw51_sla_merge.nc' ;
    tmpnc       = 'etrack_usw51_sla_tmp.nc' ;
    tracknc     = 'etrack_usw51_sla_track.nc' ;

    etrack_method= 'sla'                      ; % eddy tracking method, SLA/Q 
    data_type   = 'ROMS'                      ; % altimetry or model or ROMS
    title       = 'etrack_SLA_roms_circ_ssn_03'     ; % A general title 
    tunits      = 'seconds since 1901-01-01 00:00:00' ; % this information is passed
    torigin     = '01-JAN-1901 00:00:00';            % manually as the ROMS output
    


%------------------
% END of user input
%------------------

%profile on

    if ( min_life < 1 ) 
       error ([mfilename ':varchk'], '\n\t FATAL (%s) : min_life should be >= 1. \n', mfilename)
    end


    % number of input files and time/snapshots per input files

    if ( strcmp ( directory(end), '/' ) ~= 1 )
        directory = [directory, '/'] ;
    end 

    fid     = sprintf ('%0.4d', fid_start) ;
    reffile = [ directory, model, '_', filetype, '.', fid, '.nc' ] ;
    
    ncr = netcdf ( reffile, 'read' ) ;
         otime = ncr{'ocean_time'}(:) ;
         nsnap = length( ncr('time') ) ;
    ncr = close (ncr) ;

    nfiles = 1 + (fid_end - fid_start)/nsnap ;
    totsnap = nsnap * nfiles ;

    % initialize grid/domain parameters

    fprintf ('\t NOTE : Initializing Model grid .....\n')
    grid = func_init_grid(romsgrid,reffile, 'basic') ;

    % initialize

    count.t   = 0         ; % snapshot count (for output file), cumulative
    ncisz     = 0         ; % snapshot count (for output file), cumulative
    vdiag     = {'sla' 'vort' 'spd'} ; % valid variables are {'sla' 'q' 'w' 'vort' 'spd'} ; 

    contour_vals = [contour_vals] ; % colum vector
    ACvals       = sort(contour_vals,2,'ascend') ;   
    CCvals       = sort(contour_vals,2,'descend') ;   

    etr.cont    = contour_vals ;
    etr.minrad  = min_radius ; etr.maxrad  = max_radius ;
    etr.err     = mperr_area ;
    etr.nhann   = nhann      ; etr.tracker = mfilename  ;
    etr.zval    = zvalue     ; etr.ztype   = ztype      ; 
    etr.datadir = directory  ; etr.grid    = romsgrid    ;
    etr.fids    = fid_start  ; etr.fide    = fid_end    ;
    etr.minpts  = min_points ; etr.title   = title      ;
    etr.tunits  = tunits     ; etr.torigin = torigin    ;
    etr.tcalendar= tcalendar ; etr.method  = etrack_method ;

    flip_dim = 1 ;  % instead of (time,eid), make it (eid,time) with eid=unlimited 
    [bad vnames]   = func_create_eddy_file (tmpnc, vnames, totsnap, etr, flip_dim) ;
    
    cmd = ['cp -f ' tmpnc ' ' tracknc ];
    system(cmd) ;

    for ifl = fid_start:nsnap:fid_end ;
       fid     = sprintf ('%0.4d', ifl) ;
       filein  = [ directory, model, '_', filetype, '.', fid, '.nc' ] ;
   
       ncro = netcdf ( filein, 'read' ) ;
            otime = ncro{'ocean_time'}(:) ;
            nsnap1 = length( otime ) ; %ncro('time') ) ;
       ncro = close (ncro) ;       
       for isnap = 1:nsnap1 ;   
           count.t = count.t + 1 ; % snap count

           prog = func_get_hslice (grid, filein, isnap, ztype, zvalue)  ;

           [day,month] = func_get_roms_month(otime(isnap),tinit,tcalendar);
           ssn = floor((month-1)/3) + 1 ;
           prog.month = month ;
           prog.ssn = ssn ;

           diag      = func_get_diag_sla_geo (grid, prog, clim_file, vdiag, nhann) ;
           
           emask     = [] ;
           [antitmp emask] = func_get_eddies_circ_102010 ('anti', grid, diag, ACvals, min_radius, max_radius, mperr_area, min_points, min_slaamp) ;
           cyctmp  = func_get_eddies_circ_102010 ('cyc', grid, diag, CCvals, min_radius, max_radius, mperr_area, min_points, min_slaamp) ;
           
           ncyc  = size(cyctmp.r, 2) ;
           nanti = size(antitmp.r, 2) ;
           
           count.c  = ncyc ;   % need to know last valid i-point, to add new eddies
           count.a  = nanti ;  %       in later iterations
           
           fprintf('    Snap  = %03.0d/%03.0d     Cyclones  = %03.0d     Anticyclones  = %03.0d \n', count.t, totsnap, ncyc, nanti) ;

           ncisz = func_write_eddy_fields (tmpnc, isnap, count, cyctmp, antitmp, bad, filein, ncisz) ;

        end
     end   


