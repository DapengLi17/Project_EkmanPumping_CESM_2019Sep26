%% === readme ===

% descrip: matlab scripts plot (1) wind stress amp and curl, 
% (2)linear, nonlinear Ekman W vel, and total W vel averaged over
% winter and summer season

% update history:
% v1.0 DL 2020Feb09

% extra notes:
% ===============


%% ====set up environments====
clear all;close all;clc;

   date_str='2020Mar21';
   
infile1 = '../data_after_manipulation/TauKesaiRoWvelTempGlobalAprSepAv_2020Mar17.mat';
infile2 = '../data_after_manipulation/TauKesaiRoWvelTempGlobalOctMarAv_2020Mar17.mat';

% Constants4CESM_Global
addpath(genpath('Func4EkmanProject/'))

% % ### Kuroshio extension (KE) ###
% % lon, lat limits for contour plots
%   lat_limits = [25 55]; 
%   lat_ticks = [30:10:50];
%   lon_limits = [130 190];
%   lon_ticks = [130:20:190];
% 
% % lon, lat for text in fig  
%   Lon4text = lon_limits(1)+1.5;
%   Lat4text = lat_limits(1)+2; 
%   
% pic1 = ['../pics/WvelContoursApr2Sep_Oct2MarAvKE_',date_str,'.png'];
% % ################################

% % ### Gulf Stream (GS) ###
% % lon, lat limits for contour plots
%   lat_limits = [27 45]; 
%   lat_ticks = [30:5:45];
%   lon_limits = [275 340];
%   lon_ticks = [280:20:340];
%   
% % lon, lat for text in fig  
%   Lon4text = lon_limits(1)+1.5;
%   Lat4text = lat_limits(1)+2; 
%   
% pic1 = ['../pics/WvelContoursApr2Sep_Oct2MarAvGS_',date_str,'.png'];
% % ########################

% % ### Agulhas Return currents (AR) ###
% % lon, lat limits for contour plots
%   lat_limits = [-50 -27]; 
%   lat_ticks = [-45:5:-30];
%   lon_limits = [0 100];
%   lon_ticks = [0:20:100];
%   
% % lon, lat for text in fig  
%   Lon4text = lon_limits(1)+1.5;
%   Lat4text = lat_limits(1)+2; 
%   
% pic1 = ['../pics/WvelContoursApr2Sep_Oct2MarAvAR_',date_str,'.png'];
% % ########################

% ### Brazil Current region (BC) ###
% lon, lat limits for contour plots
  lat_limits = [-55 -20]; 
  lat_ticks = [-50:10:-20];
  lon_limits = [270 360];
  lon_ticks = [270:30:360];
  
% lon, lat for text in fig  
  Lon4text = lon_limits(1)+1.5;
  Lat4text = lat_limits(1)+2; 
  
pic1 = ['../pics/WvelContoursApr2Sep_Oct2MarAvBC_',date_str,'.png'];
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

W_ce_Oct2Mar = Oct2Mar.w_ce_TimeAv(indxLat,indxLon);
W_ce_Apr2Sep = Apr2Sep.w_ce_TimeAv(indxLat,indxLon);
W_te_Oct2Mar = Oct2Mar.w_te_TimeAv(indxLat,indxLon);
W_te_Apr2Sep = Apr2Sep.w_te_TimeAv(indxLat,indxLon);
W_le_Oct2Mar = Oct2Mar.w_le_TimeAv(indxLat,indxLon);
W_le_Apr2Sep = Apr2Sep.w_le_TimeAv(indxLat,indxLon);
W_ne_Oct2Mar = Oct2Mar.w_ne_TimeAv(indxLat,indxLon);
W_ne_Apr2Sep = Apr2Sep.w_ne_TimeAv(indxLat,indxLon);

clear Apr2Sep Oct2Mar
% --- NaN the Equatorial Regions ---
% indx4NaN = find(abs(lat_1d)<=10);
% 
% % size W_te_Oct2Mar: 1400x3600, lon_1d: 1x3600, lat_1d: 1400x1
% W_te_Apr2Sep(indx4NaN,:)=NaN; W_te_Oct2Mar(indx4NaN,:)=NaN;
% W_le_Apr2Sep(indx4NaN,:)=NaN; W_le_Oct2Mar(indx4NaN,:)=NaN;
% W_ne_Apr2Sep(indx4NaN,:)=NaN; W_ne_Oct2Mar(indx4NaN,:)=NaN;
% ======================


%% === make pics ===
f1=figure('Renderer','painters','visible','off'); % Wvel

 pic1_size = [8.5 6.5]; 
 set(f1,'units','inches','position',[5,5,pic1_size])
 cbar_limits_Wek = [-0.5 0.5];
 cbar_ticks_Wek = [-0.5:0.25:0.5];
 cbar_limits_Wce = [-1.5 1.5];
 cbar_ticks_Wce = [-1.5:0.75:1.5];
 ylabel_position=[-0.1, 0.5, 0];
 font_size = 8;
  
 % ~~~ generate subplot position ~~~  
  row_num=4;        col_num=2;
  margin_left=0.08; margin_right=0.08;
  margin_top=0.04;  margin_botm=0.06;
  pics_dist_x=0.06; pics_dist_y=0.05;
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

m_proj('miller','lon',lon_limits,'lat',lat_limits);

  subplot('Position',sbplt_posit(1,:))
   m_pcolor(lon_1d,lat_1d,W_ce_Apr2Sep.*86400);shading interp;
    polarmap;caxis(cbar_limits_Wce);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
   m_text(Lon4text,Lat4text,'a','fontsize',font_size,'FontWeight','bold');
    title('Apr-Sep');
   ylabel('W_{total}^{50m}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
   
  subplot('Position',sbplt_posit(2,:))
   m_pcolor(lon_1d,lat_1d,W_ce_Oct2Mar.*86400);shading interp;
    caxis(cbar_limits_Wce);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    title('Oct-Mar');
   m_text(Lon4text,Lat4text,'b','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(2,2) .01 sbplt_posit(2,4)-0.02], ...
       'ylim',cbar_limits_Wce,'YTick',cbar_ticks_Wce);
   title(h1,'[m/day]','FontSize',font_size-2);  
   
  subplot('Position',sbplt_posit(3,:))
   m_pcolor(lon_1d,lat_1d,W_te_Apr2Sep.*86400);shading interp;
    caxis(cbar_limits_Wek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('W_{TE}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
   m_text(Lon4text,Lat4text,'c','fontsize',font_size,'FontWeight','bold');
    
  subplot('Position',sbplt_posit(4,:))
   m_pcolor(lon_1d,lat_1d,W_te_Oct2Mar.*86400);shading interp;
    caxis(cbar_limits_Wek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
   m_text(Lon4text,Lat4text,'d','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(4,2) .01 sbplt_posit(4,4)-0.02], ...
       'ylim',cbar_limits_Wek,'YTick',cbar_ticks_Wek);
   title(h1,'[m/day]','FontSize',font_size-2);  
   
   
  subplot('Position',sbplt_posit(5,:))
   m_pcolor(lon_1d,lat_1d,W_le_Apr2Sep.*86400);shading interp;
    caxis(cbar_limits_Wek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('W_{LE}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
   m_text(Lon4text,Lat4text,'e','fontsize',font_size,'FontWeight','bold');
    
  subplot('Position',sbplt_posit(6,:))
   m_pcolor(lon_1d,lat_1d,W_le_Oct2Mar.*86400);shading interp;
    caxis(cbar_limits_Wek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
   m_text(Lon4text,Lat4text,'f','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(6,2) .01 sbplt_posit(6,4)-0.02], ...
       'ylim',cbar_limits_Wek,'YTick',cbar_ticks_Wek);
   title(h1,'[m/day]','FontSize',font_size-2);  
   

  subplot('Position',sbplt_posit(7,:))   
   m_pcolor(lon_1d,lat_1d,W_ne_Apr2Sep.*86400);shading interp;
    caxis(cbar_limits_Wek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('W_{NE}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
   m_text(Lon4text,Lat4text,'g','fontsize',font_size,'FontWeight','bold');
  
  subplot('Position',sbplt_posit(8,:))
   m_pcolor(lon_1d,lat_1d,W_ne_Oct2Mar.*86400);shading interp;
    caxis(cbar_limits_Wek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
   m_text(Lon4text,Lat4text,'h','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(8,2) .01 sbplt_posit(8,4)-0.02], ...
       'ylim',cbar_limits_Wek,'YTick',cbar_ticks_Wek);
   title(h1,'[m/day]','FontSize',font_size-2);  

   set(gcf,'color','w');
% ====================

  
%% === output data ===   
export_fig(f1,pic1,'-r200','-nocrop')
% ==================== 
