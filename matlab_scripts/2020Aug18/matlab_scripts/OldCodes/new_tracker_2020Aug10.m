tic 

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
% workdir='/scratch/user/sanjiv/Matlab-code/';

% addpath((       strcat(workdir,'PSOM/')) );
% addpath(genpath(strcat(workdir,'sanjiv_tools/')) );
% addpath(genpath(strcat(workdir,'matlabtools/')) );
% addpath(genpath(strcat(workdir,'ROMS_Jaison/gen_tools/')) );
% addpath(genpath(strcat(workdir,'ROMS_Jaison/eddy_detect/')) );
cd /scratch/user/dapengli/Projects4iHESP/Project_EkmanPumping_2019Sep27/matlab_scripts
addpath(genpath('eddy_detect'));
% addpath(genpath('CodesFromSanjiv'));

% load('hotcold_cmap','hotcold');

%---------------------------------------------------------------------------------------
    omode       = 'fast' ; 
    dbmod       = 0 ;

%    contour_vals = 15 ; 
% Kuroshio: I went with -200:200 but this might be overkill as you rarely get ssh amplitudes this high.
%    contour_vals= [-200:5:200] ; % contour values of normalized (by std of SSH) SLA (cm)
% GoM
    contour_vals= [-100:5:100] ; % contour values of normalized (by std of SSH) SLA (cm)
                                  %     to search for
    min_radius  = 20 ; %10      ; % in km , required minimum radius to define an eddy (>=)
    max_radius  = 200 ; %200      ; % in km , required maximum radius to define an eddy (>=)
    min_points  = 1 ;
    mperr_area  = 50       ; % max percentage error in area, for circle fitting to Q
    nhann       = 0        ; % apply 2D hanning smoother on vort., Q/W this many times
    min_slaamp  = 1        ; % minimum SLA amplitude in cm
    dtime       = 2 ; % spacing between snapshots 
    min_life    = 1        ; % in number of snapshots, minimum life period for eddy
    leave       = 2        ; % for joining brocken tracks

    zvalue      = -0.5      ; % z-axis/surface value, according to ztype
    ztype       = 'none' ; % depth/sigma/isopyc/none/zeta/ssh

%----Directory and grid file for the corresponding simulation---------
% Hengkai 3-km Kuroshio run
%    directory='/scratch/user/hengkai.yao/WORK/Models/CRESM_new/run/kuro03_20031015_20040330_run08/'; 
%    romsgrid    = strcat(directory,'kuro03_grid_N050_wrfmask.nc');
% Yun 3-km GoM run
%    directory='/ihesp/user/liu6/GOM_9k_nature_copernicus/Orig/';
%    romsgrid    = strcat(strrep(directory,'Orig/',''),'gom03_grd_N050_coast.nc');
directory = '../data_after_manipulation/';
romsgrid = '../data_after_manipulation/POP2ROMS_mask_2020Aug05.nc';    
%---------------------------------------------------------------------

%    clim_file   = '/home/jaison/work/store/usw51_Smean.nc' ; % Sclim
     clim_file   = []; 

%    model       = 'kuro03'   ; % prefix for input files, typically model name
    model       = 'GS'   ; % prefix for input files, typically model name
    filetype    = 'avg'     ; % average/history

%    fid_start   =  15       ; % id for first file
%    fid_end     =  25 ; %2970     ;
    fid_start   =  1       ; % id for first file
    fid_end     =  2 ; %2970     ;
    fincr       =  1;

%    flist=dir([directory 'cmpr*ocn.hi*']); 
    flist=dir([directory '*Monthly*']); 
    fno = fid_start : fincr : fid_end; 

    nfiles = length(fno);


% Kuroshio
%    tinit       = [2003 1 1] ; % model start time, for finding seasons
% GoM
    tinit       = [2010 1 1] ; % model start time, for finding seasons

    tcalendar   = '360_DAY' ;                      % 
 
%    vnames      = {'r' 'x' 'y' 'i' 'j' 'rkm' 'xkm' 'dkm' 'ykm' 'dcokm' 'area' 'slaave' 'slamax' 'cont' 'perr' 'vort' 'vortmax' 'swirlmax'} ;
     vnames      = {'r' 'x' 'y' 'i' 'j' 'rkm' 'xkm' 'dkm' 'ykm' 'dcokm' 'area' 'slaave' 'slamax' 'cont' 'perr'} ;
 
%    statsdir = '/scratch/user/sanjiv/Matlab-code/RCESM/eddy_detect/output/diag/stats/nc/GoM_3km/';
    statsdir = '../data_after_manipulation/'

    fileout     = strcat('etrack_',model,'_sla.nc') ;
    filemerge   = strcat('etrack_',model,'_sla_merge.nc') ;
    tmpnc       = strcat('etrack_',model,'_sla_tmp.nc') ;
    tracknc     = strcat('etrack_',model,'_sla_track.nc') ;

    etrack_method= 'sla'                      ; % eddy tracking method, SLA/Q 
    data_type   = 'ROMS'                      ; % altimetry or model or ROMS
    title       = 'Kuro_3km_atm_3km'     ; % A general title 
    tunits      = 'seconds since 2003-01-01 00:00:00' ; % this information is passed
    torigin     = '01-JAN-2003 00:00:00';            % manually as the ROMS output
    

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
%    reffile = [ directory, model, '_', filetype, '_', fid, '.nc' ] ;
% Kuroshio
%    reffile = [ directory,'KE03.ocn.hi.2003-12-01_06:00:00.nc' ]';
% GoM
%    reffile = [ directory,'GOM_9k_nature_copernicus.ocn.hi.2018-12-23_03:00:00.nc' ];
    reffile = ['../data_after_manipulation/POP2ROMS_Monthly_TAUXTAUYSSH_88-01.nc' ];

    ncr = netcdf ( reffile, 'read' ) ;
         otime = ncr{'ocean_time'}(:) ;
         nsnap = numel(otime) ;
    ncr = close (ncr) ;

%    nfiles = 1 + (fid_end - fid_start)/nsnap ;
    nfiles = fid_end - fid_start + 1; 
%    totsnap = nsnap * nfiles ;
%     totsnap = 365;
   totsnap = 59;
    
    % initialize grid/domain parameters

    fprintf ('\t NOTE : Initializing Model grid .....\n')
%    grid = func_init_grid(romsgrid,reffile, 'all') ;
    grid = func_init_grid(romsgrid,reffile, 'basic') ;
    grid.N=1;
    % initialize

    count.t   = 0         ; % snapshot count (for output file), cumulative
    ncisz     = 0         ; % snapshot count (for output file), cumulative
%    vdiag     = {'sla' 'vort' 'spd'} ; % valid variables are {'sla' 'q' 'w' 'vort' 'spd'} ; 
    vdiag     = {'sla'} ; % valid variables are {'sla' 'q' 'w' 'vort' 'spd'} ; 

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

% Files for eddy compositing
% The names should exactly match the name in eddy.<xxx> (in func_get_eddies_circ.m)
%    vlist_comp={'x','y','i','j','rkm','vort','slaave','sstave','sssave','sbl','utau','vtau','qnet'};
    vlist_comp={'x','y','i','j','rkm','slaave','taux','tauy'};

    cyc_stats_file  = [statsdir 'fid' sprintf('%03d',fid_start) '_' sprintf('%03d',fid_end) '_cyc.nc'];
    anti_stats_file = [statsdir 'fid' sprintf('%03d',fid_start) '_' sprintf('%03d',fid_end) '_anti.nc'];
%    cyc_stats_file  = strcat(statsdir,'cyc_eddy_stats.nc');
%    anti_stats_file = strcat(statsdir,'anti_eddy_stats.nc');

    if exist(cyc_stats_file)
      delete(cyc_stats_file)
    end
    if exist(anti_stats_file)
      delete(anti_stats_file)
    end

    func_create_ncfile(cyc_stats_file)  
    func_create_ncfile(anti_stats_file)  

%pause

% Outer FOR loop begins
    for ifile = 1:nfiles ;
%       fid     = sprintf ('%0.4d', ifl) ;
%       filein  = [ directory, model, '_', filetype, '_', fid, '.nc' ] ;
%       filein  = [ directory,'KE03.ocn.hi.2003-12-01_06:00:00.nc' ]';
    
       filein = [directory, flist(fno(ifile)).name]

       ncro = netcdf ( filein, 'read' ) ;
            otime = ncro{'ocean_time'}(:) ;
%            nsnap1 = length( otime ) ; %ncro('time') ) ;
            nsnap1 = length( otime ) ;
       ncro = close (ncro) ;       

       if ( isempty(nsnap) || nsnap < 1 )
          error ([mfilename ':timechk'], '\n\t FATAL (%s) : Cannot find valid time points in %s. \n', mfilename, reffile)
       end

%       for isnap = 3:3 ;
%        for isnap = 1:nsnap ;
        for isnap = 1:nsnap1 ;

           %isnap 
   
           count.t = count.t + 1 ; % snap count

           prog = func_get_hslice (grid, filein, isnap, ztype, zvalue)  ;

%           [year,month,day] = func_get_year(otime(isnap),[tinit,0,0,0]);  
%           [day,month] = func_get_roms_month(otime(isnap),tinit,tcalendar);
%           ssn = floor((month-1)/3) + 1 ;
%           prog.month = month ;
%           prog.ssn = ssn ;
%           ydm = str2num( [sprintf('%04d',year) sprintf('%02d',month) sprintf('%02d',day)] ) 
           ydm = datenum(otime(isnap))

           diag      = func_get_diag_sla_geo (grid, prog, clim_file, vdiag, nhann) ;

           diag.taux = prog.taux;
           diag.tauy = prog.tauy;
           
           emask     = [] ;
           [antitmp emask] = func_get_eddies_circ_102010 ('anti', grid, diag, ACvals, min_radius, max_radius, mperr_area, min_points, min_slaamp) ;
           cyctmp  = func_get_eddies_circ_102010 ('cyc', grid, diag, CCvals, min_radius, max_radius, mperr_area, min_points, min_slaamp) ;
           
%           ncyc  = size(cyctmp.r, 2) ;
%           nanti = size(antitmp.r, 2) ;
           
           ncyc  = cyctmp.n;
           nanti = antitmp.n;

           count.c  = ncyc ;   % need to know last valid i-point, to add new eddies
           count.a  = nanti ;  %       in later iterations
           
           fprintf('    Snap  = %03.0d/%03.0d     Cyclones  = %03.0d     Anticyclones  = %03.0d \n', count.t, totsnap, ncyc, nanti) ;

           if (ncyc > 0)
             func_write_eddy_stats( cyc_stats_file,cyctmp, vlist_comp,ydm,otime(isnap) ); 
           end
           if (nanti > 0)
             func_write_eddy_stats(anti_stats_file,antitmp,vlist_comp,ydm,otime(isnap) ); 
           end
%           func_write_eddy_stats( cyc_stats_file,cyctmp, vlist_comp,fid_start,otime(isnap) ); 
%           func_write_eddy_stats(anti_stats_file,antitmp,vlist_comp,fid_start,otime(isnap) ); 

%           ncisz = func_write_eddy_fields (tmpnc, isnap, count, cyctmp, antitmp, bad, filein, ncisz) ;

        end
     end   

     subplot(2,1,2);  
%      colormap(blue_red) ;
%      set (gcf,'Position', [20 470 620 650]) ;
      [hC hC] = contourf(grid.lon,grid.lat,diag.sla,contour_vals)  ; hold on ;
%      mmap_plot(grid.lon,grid.lat,diag.sla); hold on;
%      set (hC, 'LineStyle', 'none') ; axis equal ; 
%      axis ([-83 -8 -9 23]) ;
      caxis([contour_vals(1) contour_vals(end)]) ; % very important
      colorbar ;
 
      for ied = 1:antitmp.n
        [anti_cx anti_cy]  = func_get_circle (antitmp.x(ied), antitmp.y(ied), antitmp.r(ied)) ;
        %m_plot (anti_cx, anti_cy, '-k','LineWidth', 1.75)
        plot (anti_cx, anti_cy, '-k','LineWidth', 1.75)
      end
      
      for ied = 1:cyctmp.n
        [cyc_cx cyc_cy]  = func_get_circle (cyctmp.x(ied), cyctmp.y(ied), cyctmp.r(ied)) ;
        %m_plot (cyc_cx, cyc_cy, '-m','LineWidth',1.75)
        plot (cyc_cx, cyc_cy, '-m','LineWidth',1.75)
      end

     hold off 









toc  
