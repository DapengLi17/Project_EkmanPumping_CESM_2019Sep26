%
%  Description : Matlab script to find statistics of eddies from eddy
%                   tracking outputs.
%
%  Input :  ROMS grid file
%           Eddy tracking output NetCDF file
%           criterion for selecting eddies
%
%  Output : Eddy statistics like
%              - population     
%              - mean radius    (km)
%              - mean life time (days)
%              - eddy interior averages
%
%  Modifications : The original script wrote for idealized upwelling/coupled
%                     experiment has been modified for USW51 solutions.
%                     Dec/09/2008
%

%_____________________________________________________________________________
% This is a part of Eddy Tracking Tool Kit, Copyright @ UCLA ROMS Team, 2008
% Written By : Jaison Kurian (jaison@atmos.ucla.edu)
% Written On : Nov/05/2008
%_____________________________________________________________________________

%--------------------------
    close all ; % clear all ; 
    global  omode dbmod
%--------------------------

%-------------------
% USER INPUTS
%-------------------

    dbmod       = 0      ;                                   % debugging mode 

    romsgrid    = '/home/jaison/work/store/usw51_grd.nc'; % model grid
    omode       = 'fast' ;                                   % operation mode 

    reffile     = '/mnt/shirma/capet/USW5.1/CLIM_QSCATXA_SODA_HFLUX/usw51_avg.0000.nc' ;
    eddy_file   = './store_eff/etrack_usw51_sla.nc';
   
    vnames      = {'x' 'y' 'i' 'rkm' 'dcokm'} ; % variables
                         % to handle from eddy tracking output file

    min_snap    =  181; %45   ; % minimum life in number of snapshots/records 
    max_dfco    = 1E10 ; % maximum distance from the coast in km
    min_lat     =  20 ; %34  ; % minimum Longitude (deg_E); applicable for eddy init point
    max_lat     =  70 ; %46  ; % maximum Latitude (deg_N); applicable for eddy init point
    tids        =   1  ; % starting time id for period of interest
    tide        = 1529 ; %2895 ; %1530 ; % ending     "

    dtime       = 2    ; % spacing between two snapshots
    mark_birth  = true ;

    pdf         = true ; % true for pdf output and false for not


%-------------------
% END OF USER INPUTS
%-------------------

    nvar        = length(vnames) ;

    % initialize grid parameters/variables
    fprintf ('\t  NOTE : Initializing model grid ....\n')

   % grid = func_init_grid(romsgrid, reffile) ;

  % grid 
  nc = netcdf (romsgrid,'read') ;
    lon  = nc{'lon_rho'}(:) ;
    lat  = nc{'lat_rho'}(:) ;
    mask = nc{'mask_rho'}(:) ;
  nc = close (nc) ;
  mask(mask==0) = NaN ;
  [jsz isz] = size(lon) ;   

  % Model domain

   [xdom ydom ] = func_domain (lon, lat) ;

    % initialize eddy tracking file
       % tids and tide is not used here  as it will affect the 
       %   estimation of mean lifetime of eddies present in this time period.

       % this is the most time cosuming part...so do it only if required.

   if ( exist ('cyc', 'var') ~= 1 ) 
      fprintf ('\t  NOTE : Initializing CYC ....\n')
      cyc  = func_init_eddy ('cyc', eddy_file, vnames) ; 
   end
   if ( exist ('anti', 'var') ~= 1 ) 
      fprintf ('\t  NOTE : Initializing ANTI ....\n')
      anti = func_init_eddy ('anti', eddy_file, vnames) ; 
   end

    % find statistics
       % tids and tide is used here to find average stat over this period 
     fprintf ('\t  NOTE : Picking CYC ....\n')

    [cid cstat cbx cby cdx cdy cbt clife] = func_pick_eddy (cyc, min_snap, max_dfco, min_lat, max_lat, tids, tide) ;
     fprintf ('\t  NOTE : Picking ANTI ....\n')
    [aid astat abx aby adx ady abt alife] = func_pick_eddy (anti, min_snap, max_dfco, min_lat, max_lat, tids, tide) ;
    
     pc = 100 * cstat.n/(cstat.n+astat.n) ;
     pa = 100 * astat.n/(cstat.n+astat.n) ;

     fprintf ('\n\t Min_Life = %5i Days   Max_dfco = %5i km  Lat_min = %i  Lat_max = %i\n', (min_snap-1)*dtime, max_dfco, min_lat, max_lat) 
     fprintf ('\t cyc  =  %5i (%5.2f%%)  Anti =  %5i (%5.2f%%)\n', cstat.n, pc, astat.n, pa) 

     cx = cyc.x(:,cid) ;
     cy = cyc.y(:,cid) ;
     ax = anti.x(:,aid) ;
     ay = anti.y(:,aid) ;
      
     latmin = 20 ; 
     latmax = 50 ;
     lonmin = 216 ;
     lonmax = 248 ; 
     fontsize = 14 ; 

 
     fprintf ('\t  NOTE : Making figures ....\n')

     figure (11) ; set (gcf,'Position', [10 300 1000 600]) ;
     axes('position',[0.07 0.07 0.42 0.84]);
        m_proj('mercator','lon',[lonmin lonmax],'lat',[latmin latmax]);
        
        ax(end+1,:) = NaN ;
        ay(end+1,:) = NaN ; 
        m_plot(ax,ay,'r');  hold on ;
        if ( mark_birth ) 
           m_plot(abx,aby,'ok','MarkerFaceColor','k','MarkerSize',3.5) ;% hold on
        end
    %   m_plot(adx,ady,'^g','MarkerFaceColor','g','MarkerSize',3.5)
    %   m_plot(xdom,ydom,'-k','LineWidth', 0.5);  % place it here so that, no lines
                                                  %   over land
        
        %m_grid('box','fancy','xtick',[lonmin:10:lonmax],'ytick',[latmin:10:latmax],'fontsize',fontsize);
        m_grid('xtick',[lonmin:10:lonmax],'ytick',[latmin:10:latmax],'fontsize',fontsize);
        m_coast('patch',[0.9 0.9 0.9]);
        %title(['ANTI : Dist = ',num2str(max_dfco),'km, Min\_life = ',num2str(min_snap*2), ' days, Count=',num2str(astat.n),', ',(sprintf('%d',round(pa))),'%'],'fontsize',fontsize);
        title(['ANTI : Min\_life = ',num2str((min_snap-1)*dtime), ' days, Count=',num2str(astat.n),', ',(sprintf('%d',round(pa))),'%'],'fontsize',fontsize);

        hold off ;

    axes('position',[0.54 0.07 0.42 0.84]);
        
        cx(end+1,:) = NaN ;
        cy(end+1,:) = NaN ; 
        m_plot(cx,cy,'b');      hold on ; 
        if ( mark_birth ) 
           m_plot(cbx,cby,'ok','MarkerFaceColor','k','MarkerSize',3.5) ;% hold on% for filled black circles
        end
   %    m_plot(cdx,cdy,'^g','MarkerFaceColor','g','MarkerSize',3.5)
   %    m_plot(xdom,ydom,'-k','LineWidth', 0.5); 
        
        m_grid('xtick',[lonmin:10:lonmax],'ytick',[latmin:10:latmax],'fontsize',fontsize);
        m_coast('patch',[0.9 0.9 0.9]);
        %title(['CYC  : Dist = ',num2str(max_dfco),'km, Min\_life = ',num2str(min_snap*2), ' days, Count=',num2str(cstat.n),', ',(sprintf('%d',round(pc))),'%'],'fontsize',fontsize);
        title(['CYC  : Min\_life = ',num2str((min_snap-1)*dtime), ' days, Count=',num2str(cstat.n),', ',(sprintf('%d',round(pc))),'%'],'fontsize',fontsize);
    
        hold off ;
        %print -djpeg100 fig_cyc.jpg

        if ( pdf ) 
          pdfname = [mfilename, '_', num2str(max_dfco), 'km_',num2str((min_snap-1)*dtime),'days.pdf'] ;
          if (max_dfco > 1E4) 
             pdfname = [mfilename, '_', num2str((min_snap-1)*dtime),'days.pdf'] ;
          end 
          
          save2pdf (pdfname, 11, 600)
        end
