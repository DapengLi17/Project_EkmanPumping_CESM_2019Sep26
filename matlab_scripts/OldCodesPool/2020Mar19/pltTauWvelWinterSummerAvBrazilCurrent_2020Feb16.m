%% === readme ===

% descrip: matlab scripts plot (1) wind stress amp and curl, 
% (2)linear, nonlinear Ekman W vel, and total W vel averaged over
% winter and summer season

% update history:
% v1.0 DL 2020Feb09

% extra notes:
% ===============


%% === set up environments ===
clear all;close all;clc;

   date_str='2020Feb16';

infile1 = '../data_after_manipulation/TauKesaiRoWvelTempGlobalBorealSummerAv_2020Feb08.mat';
infile2 = '../data_after_manipulation/TauKesaiRoWvelTempGlobalBorealWinterAv_2020Feb08.mat';

pic1 = ['../pics/TauAmpCurlWinterSummerAvBrazilCurrent_',date_str,'.png'];
pic2 = ['../pics/WvelEkmanCESMWinterSummerAvBarzilCurrent_',date_str,'.png'];

   addpath(genpath('Func4EkmanProject/'))
  
% Constants4CESM_Global

   lat_limits = [-55 -30];  
   lon_limits = [290 320];
   lon_ticks = [290:10:320];
   lat_ticks = [-50:10:30];
%=============================


%% === load data ===
% only for southern hemishpere since its winter is boreal summer
Wt = load(infile1);
Su = load(infile2); 
%===================


%% === data analysis ===
lat_1d_raw = Su.lat_1d;
lon_1d_raw = Su.lon_1d;

indxLat = find(lat_1d_raw >= lat_limits(1) & lat_1d_raw <= lat_limits(2));
lat_1d = lat_1d_raw(indxLat);
indxLon = find(lon_1d_raw >= lon_limits(1) & lon_1d_raw <= lon_limits(2));
lon_1d = lon_1d_raw(indxLon);

tau_amp_su = Su.tau_amp_TimeAv(indxLat,indxLon);
tau_amp_wt = Wt.tau_amp_TimeAv(indxLat,indxLon);
tau_curl_su = Su.tau_curl_TimeAv(indxLat,indxLon);
tau_curl_wt = Wt.tau_curl_TimeAv(indxLat,indxLon);
W_li_wt = Wt.w_li_TimeAv(indxLat,indxLon);
W_li_su = Su.w_li_TimeAv(indxLat,indxLon);
W_nl_wt = Wt.w_nl_TimeAv(indxLat,indxLon);
W_nl_su = Su.w_nl_TimeAv(indxLat,indxLon);
W_ce_wt = Wt.w_ce_TimeAv(indxLat,indxLon);
W_ce_su = Su.w_ce_TimeAv(indxLat,indxLon);

clear Wt Su
% ======================


%% === make pics ===
f1=figure('Renderer','painters','visible','off'); % wind stress amp

  pic1_size = [8.5 4]; 
  set(f1,'units','inches','position',[5,5,pic1_size])
  cbar_limits_TauAmp  = [0 0.3];
  cbar_ticks_TauAmp  = [0:0.1:0.3];
  cbar_limits_TauCurl = [-4 4];
  cbar_ticks_TauCurl = [-4:2:4];
  x4text = -2.9;
  y4text = -1.3;
  ylabel_position=[-0.1, 0.5, 0];
  font_size = 8;
  
% ~~~ generate subplot position ~~~  
  row_num=2;        col_num=2;
  margin_left=0.08; margin_right=0.08;
  margin_top=0.05;  margin_botm=0.08;
  pics_dist_x=0.06; pics_dist_y=0.06;
 
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

m_proj('miller','lon',lon_limits,'lat',lat_limits);

 subplot('Position',sbplt_posit(1,:))
   m_pcolor(lon_1d,lat_1d,tau_amp_su);shading interp;polarmap;
    caxis(cbar_limits_TauAmp);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
% --- save codes below as reference ---
%     set(gca,'color',[0.7,0.7,0.7], ... % set ax background color as gray for NaN (land)
%       'xlim',lon_plt_limits,'xtick',lon_plt_ticks,'XTickLabel',lon_plt_tickslabel, ...
%       'ylim',lat_plt_limits,'ytick',lat_plt_ticks,'YTickLabel',lat_plt_tickslabel, ...
%       'Fontsize',font_size);
% -------------------------------------
    title('Summer Average');
    ylabel('\tau','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'a','fontsize',font_size,'FontWeight','bold');
    
 subplot('Position',sbplt_posit(2,:))
   m_pcolor(lon_1d,lat_1d,tau_amp_wt);shading interp;
    caxis(cbar_limits_TauAmp);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    title('Winter Average');
    text(x4text,y4text,'b','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(2,2) .01 sbplt_posit(2,4)-0.02], ...
       'ylim',cbar_limits_TauAmp,'YTick',cbar_ticks_TauAmp);
   title(h1,'[N/m^2]','FontSize',font_size-2);    
   
   
 subplot('Position',sbplt_posit(3,:))
   m_pcolor(lon_1d,lat_1d,tau_curl_su.*10^7);shading interp;
    caxis(cbar_limits_TauCurl);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('Curl(\tau)','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'c','fontsize',font_size,'FontWeight','bold');
    
 subplot('Position',sbplt_posit(4,:))
   m_pcolor(lon_1d,lat_1d,tau_curl_wt.*10^7);shading interp;
    caxis(cbar_limits_TauCurl);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    text(x4text,y4text,'d','fontsize',font_size,'FontWeight','bold');
% add colorbar 
    h1=colorbar('FontSize',font_size-2);
    set(h1,'position',[0.95 sbplt_posit(4,2) .01 sbplt_posit(4,4)-0.02], ...
       'ylim',cbar_limits_TauCurl,'YTick',cbar_ticks_TauCurl);
    title(h1,'[x10^{-7} N/m^3]','FontSize',font_size-3);  
% ========================


%% === make pics ===
f2=figure('Renderer','painters','visible','off'); % Wvel

 pic2_size = [8.5 6.5]; 
 set(f2,'units','inches','position',[5,5,pic2_size])
 cbar_limits_Wek = [-0.4 0.4];
 cbar_ticks_Wek = [-0.4:0.2:0.4];
 cbar_limits_Wce = [-1.2 1.2];
 cbar_ticks_Wce = [-1.2:0.6:1.2];
 x4text = -2.9;
 y4text = -1.3;
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
   m_pcolor(lon_1d,lat_1d,W_ce_su.*86400);shading interp;polarmap;
    caxis(cbar_limits_Wce);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
% --- save codes below as reference ---
%     set(gca,'color',[0.7,0.7,0.7], ... % set ax background color as gray for NaN (land)
%       'xlim',lon_plt_limits,'xtick',lon_plt_ticks,'XTickLabel',lon_plt_tickslabel, ...
%       'ylim',lat_plt_limits,'ytick',lat_plt_ticks,'YTickLabel',lat_plt_tickslabel, ...
%       'Fontsize',font_size);
% -------------------------------------
    title('Summer Average');
    ylabel('W_{eddy}^{50m}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'a','fontsize',font_size,'FontWeight','bold');
   
  subplot('Position',sbplt_posit(2,:))
   m_pcolor(lon_1d,lat_1d,W_ce_wt.*86400);shading interp;
    caxis(cbar_limits_Wce);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    title('Winter Average');
    text(x4text,y4text,'b','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(2,2) .01 sbplt_posit(2,4)-0.02], ...
       'ylim',cbar_limits_Wce,'YTick',cbar_ticks_Wce);
   title(h1,'[m/day]','FontSize',font_size-2);  
   
   
  subplot('Position',sbplt_posit(3,:))
   m_pcolor(lon_1d,lat_1d,W_nl_su.*86400);shading interp;
    caxis(cbar_limits_Wek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('W_{\xi}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'c','fontsize',font_size,'FontWeight','bold');
   
  subplot('Position',sbplt_posit(4,:))
   m_pcolor(lon_1d,lat_1d,W_nl_wt.*86400);shading interp;
    caxis(cbar_limits_Wek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    text(x4text,y4text,'d','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(4,2) .01 sbplt_posit(4,4)-0.02], ...
       'ylim',cbar_limits_Wek,'YTick',cbar_ticks_Wek);
   title(h1,'[m/day]','FontSize',font_size-2);  
   
   
  subplot('Position',sbplt_posit(5,:))
   m_pcolor(lon_1d,lat_1d,W_li_su.*86400);shading interp;
    caxis(cbar_limits_Wek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('W_{f}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'e','fontsize',font_size,'FontWeight','bold');
    
  subplot('Position',sbplt_posit(6,:))
   m_pcolor(lon_1d,lat_1d,W_li_wt.*86400);shading interp;
    caxis(cbar_limits_Wek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    text(x4text,y4text,'f','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(6,2) .01 sbplt_posit(6,4)-0.02], ...
       'ylim',cbar_limits_Wek,'YTick',cbar_ticks_Wek);
   title(h1,'[m/day]','FontSize',font_size-2);  
   

  subplot('Position',sbplt_posit(7,:))   
   m_pcolor(lon_1d,lat_1d,(W_nl_su-W_li_su).*86400);shading interp;
    caxis(cbar_limits_Wek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    ylabel('W_{\xi} - W_{f}','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'g','fontsize',font_size,'FontWeight','bold');
  
  subplot('Position',sbplt_posit(8,:))
   m_pcolor(lon_1d,lat_1d,(W_nl_wt-W_li_wt).*86400);shading interp;
    caxis(cbar_limits_Wek);
   m_coast('patch',[.7 .7 .7],'edgecolor','none');
   m_grid('linestyle','none','Fontsize',font_size, ...
       'xtick',lon_ticks,'ytick',lat_ticks);axis normal;
    text(x4text,y4text,'h','fontsize',font_size,'FontWeight','bold');
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.95 sbplt_posit(8,2) .01 sbplt_posit(8,4)-0.02], ...
       'ylim',cbar_limits_Wek,'YTick',cbar_ticks_Wek);
   title(h1,'[m/day]','FontSize',font_size-2);  
% ====================

  
%% === output data ===   
export_fig(f1,pic1,'-r200','-nocrop')
export_fig(f2,pic2,'-r200','-nocrop')
% ==================== 
