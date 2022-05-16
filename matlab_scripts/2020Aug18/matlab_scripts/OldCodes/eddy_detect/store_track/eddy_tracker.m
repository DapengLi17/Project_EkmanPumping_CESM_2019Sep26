%
%  Description : Eddy tracker, SLA-based method for ROMS.
%
%                  Please see the README file for the details of this version.
%
%  Updated : July/01/2009
%            Oct/06/2010 : see eddy_viewer for updates.
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

    romsgrid    = '/home/jaison/work/store/usw51_grd.nc'; 

    tinit       = [1901 1 1] ; % model start time, for finding seasons
    tcalendar   = '360_DAY' ;                      % 
 
    vnames      = {'r' 'x' 'y' 'i' 'j' 'rkm' 'xkm' 'dkm' 'ykm' 'dcokm' 'area' 'slaave' 'slamax' 'cont' 'perr' 'vort' 'vortmax' 'swirlmax'} ;
      
    
    fileout     = 'etrack_usw51_sla.nc' ;
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

    % initialize grid/domain parameters

    if ( strcmp ( directory(end), '/' ) ~= 1 )
        directory = [directory, '/'] ;
    end 

    fid     = sprintf ('%0.4d', fid_start) ;
    reffile = [ directory, model, '_', filetype, '.', fid, '.nc' ] ;
    fprintf ('\t NOTE : Initializing Model grid .....\n')
    grid = func_init_grid(romsgrid,reffile, 'basic') ;

    % create output NetCDF file and write header part  

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

   
    fprintf ('\t NOTE : Tracking Eddies .....\n')
    Atrack = func_track_eddy_20090825_rpr ('anti', min_life, tmpnc ) ;
    Ctrack = func_track_eddy_20090825_rpr ('cyc', min_life, tmpnc ) ;

    bad =  -1.e+34 ;
    fprintf ('\t NOTE : Writing Eddy Tracks ........\n')
    func_write_eddy_tracks (Ctrack, Atrack, tracknc, tmpnc, vnames, bad)
  
    fprintf ('\t NOTE : Applying Minimum Life........\n')
    func_min_life (grid, min_life, etr, tracknc, fileout, bad ) ; 
 
    fprintf ('\n\t DONE !!!!! ........\n')
 
