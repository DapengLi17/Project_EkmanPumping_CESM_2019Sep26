%% ===readme===

% descrip: matlab scripts plot W vel during winter storms from 
% 2) Stern nonlinear Ekman vertical pumping velocity
% 1) linear Ekman vertical pumping velocity
% 3) CESM model output

% update history:
% v1.0 DL 2019Dec05
% v1.1 DL 2019Dec09

% extra notes:
% =============


%% ====set up environments====
clear all;close all;clc;

   date_str='2020Mar21';
   
infile1 = '../data_after_manipulation/QekmanEddyGlobalAprSepAv_2020Mar18.mat';
infile2 = '../data_after_manipulation/QekmanEddyGlobalOctMarAv_2020Mar18.mat';

% Constants4CESM_Global
addpath(genpath('Func4EkmanProject/'))

% % ### Kuroshio extension (KE) ###
% % lon, lat limits for contour plots
%   lat_limits = [25 55]; 
%   lat_ticks = [30:10:50];
%   lon_limits = [130 190];
%   lon_ticks = [130:20:190];
% 
% % areas for av dpth prof and PDF  
%   bndry_lon = [140 140 170 170 140];
%   bndry_lat = [45 30 30 45 45];
% 
% % lon, lat for test in fig  
%   Lon4text = lon_limits(1)+1.5;
%   Lat4text = lat_limits(1)+2; 
%   
% pic1 = ['../pics/QekmanEddyContoursApr2Sep_Oct2MarAvKE_',date_str,'.png'];
% % ################################

% % ### Gulf Stream (GS) ###
% % lon, lat limits for contour plots
%   lat_limits = [27 45]; 
%   lat_ticks = [30:5:45];
%   lon_limits = [275 340];
%   lon_ticks = [280:20:340];
% 
% % areas for av dpth prof and PDF  
%   bndry_lon = [285 285 305 305 285];
%   bndry_lat = [42 35 35 42 42];
%   
% % lon, lat for test in fig  
%   Lon4text = lon_limits(1)+1.5;
%   Lat4text = lat_limits(1)+2; 
%   
% pic1 = ['../pics/QekmanEddyContoursApr2Sep_Oct2MarAvGS_',date_str,'.png'];
% % ########################

% % ### Agulhas Return currents (AR) ###
% % lon, lat limits for contour plots
%   lat_limits = [-50 -27]; 
%   lat_ticks = [-45:5:-30];
%   lon_limits = [0 100];
%   lon_ticks = [0:20:100];
% 
% % areas for av dpth prof and PDF  
%   bndry_lon = [10 10 60 60 10];
%   bndry_lat = [-35 -45 -45 -35 -35];
%   
% % lon, lat for test in fig  
%   Lon4text = lon_limits(1)+1.5;
%   Lat4text = lat_limits(1)+2; 
%   
% pic1 = ['../pics/QekmanEddyContoursApr2Sep_Oct2MarAvAR_',date_str,'.png'];
% % ########################

% ### Brazil Current region (BC) ###
% lon, lat limits for contour plots
  lat_limits = [-55 -20]; 
  lat_ticks = [-50:10:-20];
  lon_limits = [270 360];
  lon_ticks = [270:30:360];

% areas for av dpth prof and PDF  
  bndry_lon = [302 302 315 315 302];
  bndry_lat = [-35 -50 -50 -35 -35];
  
% lon, lat for test in fig  
  Lon4text = lon_limits(1)+1.5;
  Lat4text = lat_limits(1)+2; 
  
pic1 = ['../pics/QekmanEddyContoursApr2Sep_Oct2MarAvBC_',date_str,'.png'];
% ########################
%=============================


%% === load data ===
Apr2Sep = load(infile1);
Oct2Mar = load(infile2);
% ==================


%% === data analysis ===
lat_1d_raw = Apr2Sep.lat_1d;
lon_1d_raw = Apr2Sep.lon_1d;

indxLat = find(lat_1d_raw >= lat_limits(1) & lat_1d_raw <= lat_limits(2));
lat_1d = lat_1d_raw(indxLat);
indxLon = find(lon_1d_raw >= lon_limits(1) & lon_1d_raw <= lon_limits(2));
lon_1d = lon_1d_raw(indxLon);

Q_te_Oct2Mar = Oct2Mar.Q_te(indxLat,indxLon);
Q_te_Apr2Sep = Apr2Sep.Q_te(indxLat,indxLon);
Q_le_Oct2Mar = Oct2Mar.Q_le(indxLat,indxLon);
Q_le_Apr2Sep = Apr2Sep.Q_le(indxLat,indxLon);
Q_ne_Oct2Mar = Oct2Mar.Q_ne(indxLat,indxLon);
Q_ne_Apr2Sep = Apr2Sep.Q_ne(indxLat,indxLon);
Q_ce_Oct2Mar = Oct2Mar.Q_ce(indxLat,indxLon);
Q_ce_Apr2Sep = Apr2Sep.Q_ce(indxLat,indxLon);

clear Apr2Sep Oct2Mar

% --- NaN the Equatorial Regions ---
% indx4NaN = find(abs(lat_1d)<=10);
% 
% % size W_te_Oct2Mar: 1400x3600, lon_1d: 1x3600, lat_1d: 1400x1
% Q_te_Apr2Sep(indx4NaN,:)=NaN; Q_te_Oct2Mar(indx4NaN,:)=NaN;
% Q_le_Apr2Sep(indx4NaN,:)=NaN; Q_le_Oct2Mar(indx4NaN,:)=NaN;
% Q_ne_Apr2Sep(indx4NaN,:)=NaN; Q_ne_Oct2Mar(indx4NaN,:)=NaN;
% ======================


%% === make pics ===
f1=figure('Renderer','painters','visible','off');

  pic1_size = [8.5 6.5]; 
  set(f1,'units','inches','position',[5,5,pic1_size])
%   set(f1,'units','normalized','position',[0,0,1,1])
  font_size = 8;
  cbar_limits_Qce = [-150 150];
  cbar_ticks_Qce = [-150:75:150];
  cbar_limits_Qek = [-50 50]; 
  cbar_ticks_Qek = [-50:25:50];
  ylabel_position=[-0.1, 0.5, 0];
 
 % ~~~ generate subplot position ~~~  
  row_num=4;col_num=2;
  margin_left=0.08;margin_right=0.08;
  margin_top=0.04;margin_botm=0.06;
  pics_dist_x=0.05; pics_dist_y=0.04;
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

m_proj('miller','lon',lon_limits,'lat',lat_limits);
   
  subplot('Position',sbplt_posit(1,:))
   m_pcolor(lon_1d,lat_1d,Q_ce_Apr2Sep);shading interp;
    polarmap;caxis(cbar_limits_Qce);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    title('Apr-Sep');
    ylabel('Q_{total}^{50m}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
   m_text(Lon4text,Lat4text,'a','fontsize',font_size,'FontWeight','bold');
   m_line(bndry_lon,bndry_lat,'linewi',1.5,'color','k');
   
  subplot('Position',sbplt_posit(2,:))
   m_pcolor(lon_1d,lat_1d,Q_ce_Oct2Mar);shading interp;
    caxis(cbar_limits_Qce);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    title('Oct-Mar');
   m_text(Lon4text,Lat4text,'b','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(2,2) .01 sbplt_posit(2,4)-0.02], ...
       'ylim',cbar_limits_Qce,'YTick',cbar_ticks_Qce);
   title(h1,'[W/m^2]','FontSize',font_size-2);  
   
  subplot('Position',sbplt_posit(3,:))
   m_pcolor(lon_1d,lat_1d,Q_te_Apr2Sep);shading interp;
    caxis(cbar_limits_Qek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('Q_{TE}^{50m}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
   m_text(Lon4text,Lat4text,'c','fontsize',font_size,'FontWeight','bold');
   
  subplot('Position',sbplt_posit(4,:))
   m_pcolor(lon_1d,lat_1d,Q_te_Oct2Mar);shading interp;
    caxis(cbar_limits_Qek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
   m_text(Lon4text,Lat4text,'d','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(4,2) .01 sbplt_posit(4,4)-0.02], ...
       'ylim',cbar_limits_Qek,'YTick',cbar_ticks_Qek);
   title(h1,'[W/m^2]','FontSize',font_size-2);  

  subplot('Position',sbplt_posit(5,:))
   m_pcolor(lon_1d,lat_1d,Q_le_Apr2Sep);shading interp;
    caxis(cbar_limits_Qek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('Q_{LE}^{50m}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
   m_text(Lon4text,Lat4text,'e','fontsize',font_size,'FontWeight','bold');
   
  subplot('Position',sbplt_posit(6,:))
   m_pcolor(lon_1d,lat_1d,Q_le_Oct2Mar);shading interp;
    caxis(cbar_limits_Qek); 
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
   m_text(Lon4text,Lat4text,'f','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(6,2) .01 sbplt_posit(6,4)-0.02], ...
       'ylim',cbar_limits_Qek,'YTick',cbar_ticks_Qek);
   title(h1,'[W/m^2]','FontSize',font_size-2);     
   
  subplot('Position',sbplt_posit(7,:))
   m_pcolor(lon_1d,lat_1d,Q_ne_Apr2Sep);shading interp;
    caxis(cbar_limits_Qek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('Q_{NE}^{50m}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
   m_text(Lon4text,Lat4text,'g','fontsize',font_size,'FontWeight','bold');
      
  subplot('Position',sbplt_posit(8,:))
   m_pcolor(lon_1d,lat_1d,Q_ne_Oct2Mar);shading interp;
    caxis(cbar_limits_Qek); 
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
   m_text(Lon4text,Lat4text,'h','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(8,2) .01 sbplt_posit(8,4)-0.02], ...
       'ylim',cbar_limits_Qek,'YTick',cbar_ticks_Qek);
   title(h1,'[W/m^2]','FontSize',font_size-2);   

   set(gcf,'color','w');
% ====================

  
%% === output data ===  
export_fig(f1,pic1,'-r200','-nocrop')

% % write statistics to data file
% fid = fopen(outfile1,'w');
% 
%  fprintf(fid,['readme: generated by DL on', date_str, ...
%      ' via code pltQekmanEddyWinterSummerAvKuroshio_2020Feb14.m \n']);
%  fprintf(fid,' \n');
% 
%  fprintf(fid,['                       Summer(Apr-Sep)  Winter(Oct-Mar) \n']);
%  fprintf(fid,['                        mean+-std        mean+-std \n']);
%  
% %
%  fprintf(fid,'Q_LIEkman(50m) [W/m2]   %.2f +- %.2f     %.2f +- %.2f\n', ...
%    [nanmean(Q_li_su(:)),nanstd(Q_li_su(:)),nanmean(Q_li_wt(:)),nanstd(Q_li_wt(:))]);
%  fprintf(fid,' \n');
% 
% %
%  fprintf(fid,'Q_NLEkman(50m) [W/m2]   %.2f +- %.2f     %.2f +- %.2f\n', ...
%    [nanmean(Q_nl_su(:)),nanstd(Q_nl_su(:)),nanmean(Q_nl_wt(:)),nanstd(Q_nl_wt(:))]);
%  fprintf(fid,' \n'); 
%  
% %
%  fprintf(fid,'Q_eddy(50m)    [W/m2]   %.2f +- %.2f   %.2f +- %.2f\n', ...
%    [nanmean(Q_ce_su(:)),nanstd(Q_ce_su(:)),nanmean(Q_ce_wt(:)),nanstd(Q_ce_wt(:))]);
%  fprintf(fid,' \n');
%  
% fclose(fid);  
% ====================
