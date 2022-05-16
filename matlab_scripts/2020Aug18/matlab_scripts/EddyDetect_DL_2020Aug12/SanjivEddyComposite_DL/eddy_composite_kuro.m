tic

 clear all ; %close all 

 workdir='/scratch/user/sanjiv/Matlab-code/';

 addpath((       strcat(workdir,'PSOM/')) );
 addpath(genpath(strcat(workdir,'sanjiv_tools/')) );
 addpath(genpath(strcat(workdir,'matlabtools/')) );
 addpath(genpath(strcat(workdir,'ROMS_Jaison/gen_tools/')) );
 addpath(genpath(strcat(workdir,'ROMS_Jaison/eddy_detect/')) );

 load('hotcold_cmap','hotcold');

% directory='/scratch/user/hengkai.yao/WORK/Models/CRESM_new/run/kuro03_20031015_20040330_run08/'; 
 directory='/ihesp/user/liu6/GOM_9k_nature_copernicus/Orig/';
 
%stats_dir ='/scratch/user/sanjiv/Matlab-code/RCESM/eddy_detect/stats/';
stats_dir ='/scratch/user/sanjiv/Matlab-code/RCESM/eddy_detect/output/diag/stats/nc/GoM_3km/';

% grid_dir='/scratch/user/sanjiv/RCESM_output/kuro03_20031015_20040330_run08/fulldomain_historyfiles/';
 grid_dir='/scratch/user/sanjiv/RCESM_output/GOM_ocn3km_atm9km/';
 
% romsgrid    = strcat(directory,'kuro03_grid_N050_wrfmask.nc');

 coordfile = strcat(grid_dir,'grid_3hr.nc'); 
 grid_struc=get_gridpar(coordfile);

 ocean_time = grid_struc.ocean_time;
 ntimes=max(size(ocean_time)); 


%----READ IN GRID PARAMETERS---------
 Cs_r    = grid_struc.Cs_r;
 theta_s = grid_struc.theta_s;
 theta_b = grid_struc.theta_b;
 hc      = grid_struc.hc;
 vtrans  = grid_struc.vtrans;
 vstret  = grid_struc.vstret;

 s_rho   = grid_struc.s_rho;
 s_w     = grid_struc.s_w;  
 nz_rho  = grid_struc.nz_rho;
 nz_w    = grid_struc.nz_w; 

 bathy = grid_struc.bathy;
 pm    = grid_struc.pm;
 pn    = grid_struc.pn;
 lonr  = grid_struc.lonr;
 latr  = grid_struc.latr;
 maskr = grid_struc.maskr; 

 maskr(maskr~=1) = NaN ;

 nx = size(pm,1); ny = size(pm,2); nxy = nx*ny;  
%---------------------------------

 filetype    = 'avg' ; % average/history
 fid_start   =  15   ; % id for first file
 fid_end     =  30  ; 

% testing purposes
% fid_start   =  67       ; % id for first file
% fid_end     =  75 ; %2970     ;

 fincr       =  1;

% flist = fid_start : fincr : fid_end; 
 flist=dir(strcat(directory,'cmpr*ocn.hi*')); 
 fno = fid_start : fincr : fid_end; 

 nfiles = length(fno);

 reffile = [ directory,'KE03.ocn.hi.2003-12-01_06:00:00.nc' ];
 ncr = netcdf ( reffile, 'read' ) ;
 ot  = ncr{'ocean_time'}(:) ;
 nsnap = numel(ot) ;
 ncr = close (ncr) ;

 cyc_file = [stats_dir,'fid15_30_cyc.nc'];

 cyc_fid = ncread( cyc_file,'fid' )'; 
 cyc_tid = ncread( cyc_file,'tid' )'; 

% Initialize eddy count
 composite.ok  = 0;
 composite.nok = 0;
 composite.tot = 0;

% Multiplier of eddy radius that sets the extent of the domain where fields 
% are saved
% reddy_fac = 6; % Kuroshio
 reddy_fac = 5; % GoM 

% Vector used for permuting 2D dimension order (to call u2rho3d, etc.)
 dimorder = [2 1];

%---DOMAIN SUBSET FOR SELECTING EDDIES------
% Kuroshio domain subset for eddy statistics. Ignore eddies outside this region (to avoid land, shallow depths)
  lon_left = 144;  
  ileft=find(lonr(:,100)-lon_left > 0 ,1,'first'); iright = nx;
% ileft   = 1 ; iright = nx;

%  lat_top = 42;
%  jtop =find(latr(100,:)-lat_top <0,1,'last'); jbottom = 1; 
  jbottom = 1 ; jtop =   ny; 
%%


%--------------------------------------------

% Final scaled grid to which values are interpolated
  [xgrid,ygrid] = meshgrid( -(reddy_fac-1):0.1:(reddy_fac-1), -(reddy_fac-1):0.1:(reddy_fac-1) );

% Initialize fields to be composited 
  composite.sst = zeros(size(xgrid));  
  composite.ssh = composite.sst;  
  composite.tau = composite.sst;  
  composite.taucurl = composite.sst;  
  composite.taudiv  = composite.sst;  
  composite.dsstdx  = composite.sst;
  composite.dsstdy  = composite.sst;
  composite.gradsst = composite.sst;
  composite.qnet = composite.sst;
  composite.qlat = composite.sst;
  composite.qsen = composite.sst;
  composite.hsbl = composite.sst;


 for ifile = 1:nfiles;

% I screwed up in missing files Dec. 21, Dec. 22 and Dec. 23 while writing out the eddy statistics file. 
     if (fno(ifile)==68 |fno(ifile)==69 |fno(ifile)==70 ) 
       continue
     end 

     fname = flist(fno(ifile)).name
     
%   fid = sprintf('%0.4d',fno(ifile)); 

%   cyc_fid_start = find(cyc_fid==flist(ifile),1,'first');
%   cyc_fid_end   = find(cyc_fid==flist(ifile),1,'last');

   for isnap = 1:nsnap;

     %isnap  
     otime = ncread([directory fname],'ocean_time',isnap,1,1);
%     otime = ocean_time( (fno(ifile)-1)*nsnap + isnap );

     cyc_start = find(cyc_tid==otime,1,'first');
     if (isempty(cyc_start)==1)
       continue
     end 
     cyc_end   = find(cyc_tid==otime,1,'last');
     cyc_neddy = cyc_end - cyc_start + 1;


     ssh = squeeze(ncread( strcat(directory,fname),'zeta',[1 1 isnap],[nx ny 1],[1 1 1] ));            ssh = double(ssh);  
     sst = squeeze(ncread( strcat(directory,fname),'temp',[1 1 nz_rho isnap],[nx ny 1 1],[1 1 1 1] )); sst = double(sst);  

     qnet = squeeze(ncread( strcat(directory,fname),'shflux'  ,[1 1 isnap],[nx ny 1],[1 1 1] )); qnet = double(qnet);  
     qlat = squeeze(ncread( strcat(directory,fname),'latent'  ,[1 1 isnap],[nx ny 1],[1 1 1] )); qlat = double(qlat);  
     qsen = squeeze(ncread( strcat(directory,fname),'sensible',[1 1 isnap],[nx ny 1],[1 1 1] )); qsen = double(qsen);  
     hsbl = squeeze(ncread( strcat(directory,fname),'Hsbl'    ,[1 1 isnap],[nx ny 1],[1 1 1] )); hsbl = double(hsbl);  

     utau = squeeze(ncread( strcat(directory,fname),'sustr',[1 1 isnap],[nx-1 ny 1],[1 1 1] ));
     vtau = squeeze(ncread( strcat(directory,fname),'svstr',[1 1 isnap],[nx ny-1 1],[1 1 1] ));  
     utau = permute(u2rho_2d( permute(utau,dimorder) ), dimorder); 
     vtau = permute(v2rho_2d( permute(vtau,dimorder) ), dimorder); 

     sst  =  sst.*maskr ;
     ssh  =  ssh.*maskr ;
     utau = utau.*maskr;
     vtau = vtau.*maskr;  

     qnet = qnet.*maskr ;
     qlat = qlat.*maskr ;
     qsen = qsen.*maskr ;
     hsbl = hsbl.*maskr ;

     tau = (utau.^2+vtau.^2).^0.5;

     taucurl=get_windcurl(rho2u_2d(permute(utau,dimorder)),rho2v_2d(permute(vtau,dimorder)),pm,pn,maskr);

     dutau_dx = (utau(2:end,:)-utau(1:end-1,:)).*...
          0.5.*(pm(2:end,:)+pm(1:end-1,:));
     dvtau_dy = (vtau(:,2:end)-vtau(:,1:end-1)).*...
          0.5.*(pn(:,2:end)+pn(:,1:end-1));
     dutau_dx = permute( u2rho_2d(permute(dutau_dx,dimorder)), dimorder);
     dvtau_dy = permute( v2rho_2d(permute(dvtau_dy,dimorder)), dimorder);
     taudiv = dutau_dx + dvtau_dy; 

     dsst_dx = (sst(2:end,:)-sst(1:end-1,:)).*...
          0.5.*(pm(2:end,:)+pm(1:end-1,:));
     dsst_dy = (sst(:,2:end)-sst(:,1:end-1)).*...
          0.5.*(pn(:,2:end)+pn(:,1:end-1));
     dsst_dx = permute( u2rho_2d(permute(dsst_dx,dimorder)), dimorder);
     dsst_dy = permute( v2rho_2d(permute(dsst_dy,dimorder)), dimorder);
     grad_sst = (dsst_dx.^2 + dsst_dy.^2).^0.5; 
%
     cyc_lon = ncread(cyc_file,'x'  ,cyc_start,cyc_neddy,1);
     cyc_lat = ncread(cyc_file,'y'  ,cyc_start,cyc_neddy,1);
     cyc_i   = ncread(cyc_file,'i'  ,cyc_start,cyc_neddy,1);
     cyc_j   = ncread(cyc_file,'j'  ,cyc_start,cyc_neddy,1);
     cyc_rkm = ncread(cyc_file,'rkm',cyc_start,cyc_neddy,1);

%     cyc_sst = ncread(cyc_file,'sstave',cyc_start,cyc_neddy,1);
%     cyc_sss = ncread(cyc_file,'sssave',cyc_start,cyc_neddy,1);

%     cyc_utau = ncread(cyc_file,'utau',cyc_start,cyc_neddy,1);
%     cyc_vtau = ncread(cyc_file,'vtau',cyc_start,cyc_neddy,1);
%     cyc_qnet = ncread(cyc_file,'qnet',cyc_start,cyc_neddy,1);

% Update the total number of eddies with each snapshot
     composite.tot = composite.tot + cyc_neddy; 

% LOOP THROUGH THE EDDIES AT EACH TIMESTEP
     for ieddy=1:cyc_neddy

       eddy_rad = cyc_rkm(ieddy);
       max_rad = reddy_fac*eddy_rad; 

% Use (i,j) on ROMS grid nearest to eddy center
       i0=cyc_i(ieddy) ; j0=cyc_j(ieddy); 
%       lat0 = latr(i0,j0); lon0 = lonr(i0,j0); 
% Use actual (lon,lat) pair describing eddy center (not necessarily on the grid anymore)
%       lat0 = latr(i0,j0); lon0 = lonr(i0,j0); 
% Use (i0,j0) to locate eddy center.
       lat0 = cyc_lat(ieddy) ; lon0 = cyc_lon(ieddy); 

%        parfor l=1:nxy;
%           [dist(l) theta(l)]=sw_dist([latr(l) cyc_lat(ieddy)],[lonr(l) cyc_lon(ieddy)],'km');
%        end

%        parfor j=1:ny;
%          for i=1:nx;
%            [dist(i,j) theta(i,j)]=sw_dist([latr(i,j) latr(cyc_i(ieddy), cyc_j(ieddy))],[lonr(i,j) lonr(cyc_i(ieddy),cyc_j(ieddy))],'km');
%            [dist(i,j) theta(i,j)]=sw_dist([latr(i,j) cyc_lat(ieddy)],[lonr(i,j) cyc_lon(ieddy)],'km');
%          end
%        end 
%        dist = dist /cyc_rkm(ieddy);


        if (abs(sw_dist([lat0 latr(1    ,jbottom)],[lon0 lon0],'km')) < max_rad | ...
            abs(sw_dist([lat0 latr(1    ,jtop   )],[lon0 lon0],'km')) < max_rad | ...
            abs(sw_dist([lat0 lat0               ],[lon0 lonr(ileft, 1)],'km')) < max_rad | ...
            abs(sw_dist([lat0 lat0               ],[lon0 lonr(iright,1)],'km')) < max_rad )
%          strcat('ieddy ', sprintf('%0.4d',ieddy), ' too close to boundary!')
          continue
        end 
%
        for i=i0:nx
           dist = sw_dist([latr(i,j0) lat0],[lonr(i,j0) lon0],'km');
           if (dist/eddy_rad > reddy_fac)
             break
           end 
        end
        imax = i-1; imin = i0-(imax-i0);
        for j=j0:ny
           dist = sw_dist([latr(i0,j) lat0],[lonr(i0,j) lon0],'km');
           if (dist/eddy_rad > reddy_fac)
             break
           end 
        end
        jmax = j-1; jmin = j0-(jmax-j0);

% Discard eddies such that a point reddy_fac*R from its center falls on land or is outside the domain

        if ( imin <= ileft | imax >= iright | jmin <= jbottom | jmax >= jtop | ...
             ~isempty( find(isnan(maskr(imin:imax,jmin:jmax))==1) ) )
          continue
        end 

        xl_subset = sw_dist([latr(imin,jmin) latr(imax,jmin)],[lonr(imin,jmin) lonr(imax,jmin)],'km'); 
        yl_subset = sw_dist([latr(imin,jmin) latr(imin,jmax)],[lonr(imin,jmin) lonr(imin,jmax)],'km'); 

% Get means of variables over [imin:imax,jmin:jmax]
        utau_mean=nanmean(nanmean(utau(imin:imax,jmin:jmax))); 
        vtau_mean=nanmean(nanmean(vtau(imin:imax,jmin:jmax))); 

%
        dist=zeros(imax-imin+1,jmax-jmin+1); arg = dist; 
        parfor j = 1:jmax-jmin+1
          for i = 1:imax-imin+1
             imove = i+imin-1; jmove = j+jmin-1;  
             [dist(i,j), arg(i,j)] = sw_dist([lat0 latr(imove,jmove)],[lon0 lonr(imove,jmove)],'km');
          end
        end
        dist = dist/eddy_rad; 
        arg(arg<0) = arg(arg<0) + 360; arg = arg*pi/180;  

        xo = dist.*cos(arg); yo = dist.*sin(arg); 

% ROTATE EDDY w.r.t winds such that wind is always from west to east
       tau_arg = atan2(vtau_mean,utau_mean); if (tau_arg<0); tau_arg = tau_arg + 2*pi; end  
       arg = arg - tau_arg;
       arg(arg<0) = arg(arg<0) + 2*pi; 

        xrot = dist.*cos(arg); yrot = dist.*sin(arg); 

%        err_center(ieddy) = sw_dist([cyc_lat(ieddy) latr(cyc_i(ieddy), cyc_j(ieddy))],[cyc_lon(ieddy) lonr(cyc_i(ieddy),cyc_j(ieddy))],'km');

       [X,Y]=meshgrid(size(xrot,1),size(xrot,2));

        for i=1:size(xrot,1);
         for j=1:size(xrot,2);
            X(j,i) = xrot(i,j);
            Y(j,i) = yrot(i,j);
         end
        end

% Rotate dsst_dx and dsst_dy into along-wind and cross-wind components
        dsst_along =  dsst_dx.*cos(tau_arg) + dsst_dy.*sin(tau_arg);
        dsst_cross = -dsst_dx.*sin(tau_arg) + dsst_dy.*cos(tau_arg);



% xgrid, ygrid are the new locations (scaled by eddy radius) where I want the interpolated values
% The interpolation step
       sst_grid = griddata(X,Y,sst(imin:imax,jmin:jmax)',xgrid,ygrid,'cubic');
       ssh_grid = griddata(X,Y,ssh(imin:imax,jmin:jmax)',xgrid,ygrid,'cubic');
       tau_grid = griddata(X,Y,tau(imin:imax,jmin:jmax)',xgrid,ygrid,'cubic');
       taucurl_grid = griddata(X,Y,  taucurl(imin:imax,jmin:jmax)',xgrid,ygrid,'cubic');
       taudiv_grid  = griddata(X,Y,   taudiv(imin:imax,jmin:jmax)',xgrid,ygrid,'cubic');
       dsstdx_grid  = griddata(X,Y,  dsst_along(imin:imax,jmin:jmax)',xgrid,ygrid,'cubic');
       dsstdy_grid  = griddata(X,Y,  dsst_cross(imin:imax,jmin:jmax)',xgrid,ygrid,'cubic');
       gradsst_grid = griddata(X,Y, grad_sst(imin:imax,jmin:jmax)',xgrid,ygrid,'cubic');

       qnet_grid = griddata(X,Y,qnet(imin:imax,jmin:jmax)',xgrid,ygrid,'cubic');
       qlat_grid = griddata(X,Y,qlat(imin:imax,jmin:jmax)',xgrid,ygrid,'cubic');
       qsen_grid = griddata(X,Y,qsen(imin:imax,jmin:jmax)',xgrid,ygrid,'cubic');
       hsbl_grid = griddata(X,Y,hsbl(imin:imax,jmin:jmax)',xgrid,ygrid,'cubic');

       composite.ssh = composite.ssh + ssh_grid; 
       composite.sst = composite.sst + sst_grid; 
       composite.tau = composite.tau + tau_grid; 
       composite.taucurl = composite.taucurl + taucurl_grid; 
       composite.taudiv  = composite.taudiv  + taudiv_grid; 
       composite.dsstdx  = composite.dsstdx  + dsstdx_grid; % x = along-wind 
       composite.dsstdy  = composite.dsstdy  + dsstdy_grid; % y = cross-wind 
       composite.gradsst = composite.gradsst + gradsst_grid;

       composite.qnet = composite.qnet + qnet_grid; 
       composite.qlat = composite.qlat + qlat_grid; 
       composite.qsen = composite.qsen + qsen_grid; 
       composite.hsbl = composite.hsbl + hsbl_grid; 

 
% Update the number of eddies that meet the specified criteria 
       composite.ok = composite.ok + 1;


% figure(2);
% subplot(2,4,max(1,mod(count,9)))
% pcolor(xrot,yrot,sst(imin:imax,jmin:jmax));shading interp
% cmocean('thermal'); colorbar
% xlim([-(reddy_fac-1) reddy_fac-1]); ylim([-(reddy_fac-1) reddy_fac-1]);
% set(gca,'xtick',[-(reddy_fac-1):reddy_fac-1]); set(gca,'ytick',[-(reddy_fac-1):reddy_fac-1]);
% pbaspect([xl_subset yl_subset 1])
%
% figure(3);
% subplot(2,4,max(1,mod(count,9)))
% pcolor(xrot,yrot,tau(imin:imax,jmin:jmax));shading interp
% cmocean('thermal'); colorbar
% caxis([0 0.1]);
% xlim([-(reddy_fac-1) reddy_fac-1]); ylim([-(reddy_fac-1) reddy_fac-1]);
% set(gca,'xtick',[-(reddy_fac-1):reddy_fac-1]); set(gca,'ytick',[-(reddy_fac-1):reddy_fac-1]);
% pbaspect([xl_subset yl_subset 1])

% pause


     end   % (for ieddy) 

  end  % (for isnap)

 end  % (for ifile) 


composite.sst = composite.sst/composite.ok; 
composite.ssh = composite.ssh/composite.ok; 
composite.tau = composite.tau/composite.ok; 
composite.taucurl = composite.taucurl/composite.ok; 
composite.taudiv  = composite.taudiv/composite.ok; 
composite.dsstdx  = composite.dsstdx/composite.ok; 
composite.dsstdy  = composite.dsstdy/composite.ok; 
composite.gradsst = composite.gradsst/composite.ok; 
composite.qnet = composite.qnet/composite.ok; 
composite.qlat = composite.qlat/composite.ok; 
composite.qsen = composite.qsen/composite.ok; 
composite.hsbl = composite.hsbl/composite.ok; 

composite.xgrid = xgrid;
composite.ygrid = ygrid; 

% SAVE COMPOSITE FIELDS
outdir=['/scratch/user/sanjiv/Matlab-code/RCESM/eddy_detect/output/diag/stats/mat/GoM_3km/'];

save([outdir 'composite_fid' sprintf('%03d',fid_start) '_' sprintf('%03d',fid_end) '.mat'], 'composite'); 

% PLOTS
make_plots=0;

if (make_plots==1)

figure(2);
pcolor(xgrid,ygrid,sst_composite);shading flat;
pbaspect([1 1 1]);
colorbar;
xlim([-4 4]); ylim([-4 4]);
set(gca,'xtick',[-5:5]);
set(gca,'ytick',[-5:5]);

figure(3);
pcolor(xgrid,ygrid,tau_composite);shading flat;
pbaspect([1 1 1]);
colorbar;
xlim([-4 4]); ylim([-4 4]);
set(gca,'xtick',[-5:5]);
set(gca,'ytick',[-5:5]);

figure(4);
pcolor(xgrid,ygrid,taucurl_composite);shading flat;
pbaspect([1 1 1]);
cmocean('balance')
colorbar;
xlim([-4 4]); ylim([-4 4]);
set(gca,'xtick',[-5:5]);
set(gca,'ytick',[-5:5]);

end 


toc
