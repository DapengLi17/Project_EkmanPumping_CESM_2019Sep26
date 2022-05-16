%
%  eddy_viewer_SSH_model.m
%
%  Description : Script to visualize eddies for a given snapshot in SSH
%                from ROMS model output.
%               
%  Use         : Use this script to make choices for various fields 
%                (like contour_vals, min_radius, max_radius, min_points,
%                mperr_area, nhann, min_slaamp, subdomain, EBndry_lon etc)
%                for the use in eddy_identifier_SSH_model.m
%
%  Assumptions : - Model SSH and Mean (seasonal/annual) SSH fields are available.
%                     SLA is computed as SSH-mean_SSH.
%                - Single or multi input files (from ROMS run).
%                - ROMS SSH (variable zeta) units is m. 
%

%____________________________________________________________
%  Written By     : Jaison Kurian (jaisonk@tamu.edu)
%                   Oceanography Dept., TAMU.
%  Date           : Jan/11/2011
%  eTRACK Version : 2.00 
%  Ref. Version   : eTRACK 1.00
%  Copyright      : See etrack_manual.pdf
%  Modifications  : 08/Jul/2014: Added options for Rutgers ROMS.
%____________________________________________________________
%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    clear all ; 
%    close all  
%---------------------------------------------------------------------------------------
    global dbmod vnames sla_t
%---------------------------------------------------------------------------------------

addpath(genpath('/scratch/user/sanjiv/Matlab-code/ROMS_Jaison/'));
addpath(genpath('/scratch/user/sanjiv/Matlab-code/matlabtools/'));

load('hotcold_cmap','hotcold');
%-------------
% USER INPUT
%-------------

    pdf         = true        ; % need pdf output
    fno         = 2           ; % figure number, to save pdf output

    dbmod       = 0            ; % debug mode: 
    
    contour_vals= [-200:5:200]   ; % contour values of normalized (by std of SSH) SLA (cm)
                                 %     to search for
    min_radius  = 40           ; % in km , required minimum radius to define an eddy (>=)
    max_radius  = 150          ; % in km , required maximum radius to define an eddy (>=)
    min_life    = 1            ; % in number of snapshots, minimum life period for eddy
    min_points  = 1            ;
    mperr_area  = 50           ; % max percentage error in area, for circle fitting to Q
    nhann       = 0            ; % apply 2D hanning smoother on vort., Q/W this many times
    min_slaamp  = 1            ; % minimum SLA amplitude in cm

    dtime       = 1            ; % spacing between snapshots 

    zvalue      = 0            ; % z-axis/surface value, according to ztype
    ztype       = 'none'       ; % depth/sigma/isopyc/none/zeta/ssh/none

    vnames      = {'r' 'x' 'y' 'i' 'j' 'rkm' 'xkm' 'dkm' 'ykm' 'dcokm' 'vort' 'vortmax' 'swirlmax' } ;

    directory='/scratch/user/hengkai.yao/WORK/Models/CRESM_new/run/kuro03_20031015_20040330_run08/'; 

    romsgrid    = strcat(directory,'kuro03_grid_N050_wrfmask.nc');
%    romsgrid    = '/ihesp/sanjiv/cresm-1.0.0/test/run/Xiao_Amazon_Ric_0.3_lmdcv_1.25/gom03_grd_N050_md15m.nc';
%    romsgrid    = '/work4/maxiaohui/crcm/grid/grid_cplr9_100313.nc' ; % ROMS Grid file (for model domain)
    romsbranch  = 'Rutgers'    ; % 'Rutgers' (myroms) or 'UCLA' 
    omode       = 'basic'      ; % 'basic' or 'all' --> selection of grid variables 
    coast_ref   = 'west'       ; % ref coastline side ('east' 'west' 'south' 'north') to find
                                 %       "distance from the coast"
%    directory   = '/normal1/maxiaohui/crcm/crcm_out/cpl_ctrl_2002/roms' ;
    clim_file   = []           ; % "Annual Mean ZETA", if not available, leave empty as [].

    model       = 'roms'      ; % prefix for input files, typically model name
    filetype    = 'var3D'        ; % average/history
    fid_start   =  8710        ; % id for first file
    isnap       =   1          ; % which snapshot to read from reference file

    tinit       = [2003 1 1]   ; % model start time, for finding seasons
    tcalendar   = '360_DAY'    ; % model calendar, for finding seasons

%------------------
% END of user input
%------------------
%profile on

%------------
% Initialize
%------------

    %----------Min_life check-------------------------------------------% 
    if ( min_life < 1 ) 
       error ([mfilename ':varchk'], '\n\t FATAL (%s) : min_life should be >= 1. \n', mfilename)
    end

    %----------Reference filename/path----------------------------------% 
    if ( strcmp ( directory(end), '/' ) ~= 1 )
        directory = [directory, '/'] ;
    end 
    fid     = sprintf ('%0.4d', fid_start) ;
%    reffile = [ directory, model, '_', filetype, '_', fid, '.nc' ] ;
%    reffile = [directory 'TXGLO.ocn.hi.2007-08-16_21:00:00.nc']; 
    reffile = [directory 'KE03.ocn.hi.2003-12-31_06:00:00.nc']; 
     
    %----------Time Variable and dimension------------------------------% 
    tdim = 'ocean_time' ; 
    tvar = 'ocean_time' ;
    if ( strcmp(lower(romsbranch),'ucla') == 1 )
        tdim = 'time' ;
    end
    ncr = netcdf ( reffile, 'read' ) ;
         otime = ncr{tvar}(:) ;
         nsnap = length( ncr(tdim) ) ;
    ncr = close (ncr) ;

    nfiles  = isnap ;
    totsnap = isnap ;

    %----------Init ROMS grid-------------------------------------------% 
    fprintf ('\t NOTE : Initializing Model grid .....\n')
    grid = func_init_grid(romsgrid,reffile, omode) ;

    %----------Init Model's Prognostic vars-----------------------------% 
    fprintf ('\t NOTE : Initializing Model Varialbes ......\n')
    prog = func_get_hslice (grid, reffile, isnap, ztype, zvalue)  ;

    %----------Find Model snapshot's month and season-------------------% 
    [day,month] = func_get_roms_month(otime(isnap),tinit,tcalendar);
    ssn         = floor((month-1)/3) + 1 ;
    prog.month  = month ;
    prog.ssn    = ssn ;

    %----------Calculate Model Diagnostic Fields------------------------% 
    fprintf ('\t NOTE : Calculating diagnostic fields ......\n')
    vdiag     = {'sla' 'vort' 'spd'} ; % valid variables are {'sla' 'q' 'w' 'vort' 'spd'} ; 
    diag      = func_get_diag_sla_geo (grid, prog, clim_file, vdiag, nhann) ;

%----------------------
% Eddy Indentificaton
%----------------------
    %----------Identify Eddies for given snapshot-----------------------% 
    fprintf ('\t NOTE : Identifying Eddies......\n')
    contour_vals = [contour_vals] ; % column vector
    ACvals       = sort(contour_vals,2,'ascend') ;   
    CCvals       = sort(contour_vals,2,'descend') ;   
    emask        = [] ;
    [antipr emask] = func_get_eddies_circ_102010 ('anti', grid, diag, ACvals, min_radius, max_radius, mperr_area, min_points, min_slaamp) ;
    cycpr          = func_get_eddies_circ_102010 ('cyc', grid, diag, CCvals, min_radius, max_radius, mperr_area, min_points, min_slaamp) ;
    fprintf('\n    Snap  = %03.0d/%03.0d     Cyclones  = %03.0d     Anticyclones  = %03.0d \n', isnap, totsnap, cycpr.n, antipr.n ) ;    

%------------
% make plots
%------------

   figure(fno) ; colormap(blue_red) ;
%      set (gcf,'Position', [20 470 620 650]) ;
      [hC hC] = contourf(grid.lon,grid.lat,diag.sla,contour_vals)  ; hold on ;
      set (hC, 'LineStyle', 'none') ; axis equal ; 
%      axis ([120 200 23 50]) ;
%      axis ([-83 -8 -9 23]) ;
      caxis([contour_vals(1) contour_vals(end)]) ; % very important
      colorbar ;
      
      for ied = 1:antipr.n
        [anti_cx anti_cy]  = func_get_circle (antipr.x(ied), antipr.y(ied), antipr.r(ied)) ;
        plot (anti_cx, anti_cy, '-k','LineWidth', 1.75)
      end
      
      for ied = 1:cycpr.n
        [cyc_cx cyc_cy]  = func_get_circle (cycpr.x(ied), cycpr.y(ied), cycpr.r(ied)) ;
        plot (cyc_cx, cyc_cy, '-r','LineWidth',1.75)
      end
      
      hold off;

    %----------Save output if required----------------------------------% 
      if ( pdf ) 
         epsname = ['ed_sla.eps'] ;
         cmd1 = ['rm -f ' epsname] ;
         system (cmd1) ;
         print(fno,'-depsc2',epsname) ;
         cmd = ['epstopdf ' epsname];
         system(cmd) ;
         system (cmd1) ;
      end

%-------
% DONE
%-------
