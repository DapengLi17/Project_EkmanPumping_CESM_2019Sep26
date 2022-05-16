%% ===readme===

% descrip: matlab scripts plot Nonlinear Ekman vertical velocity (W)
% based on Stern (1965) formulation  
% 2) Stern nonlinear Ekman vertical pumping velocity
% 1) linear Ekman vertical pumping velocity
% 3) CESM model output

% update history:
% v1.0 DL 2019Dec05
% v1.1 DL 2019Dec09
% v1.2 DL 2020Feb05

% extra notes:
% =============


%% === set up environments ===
clear all;close all;clc;

  date_str='2020Feb05'; 
  
infile1 = '../raw_data/CESM_GlobalRegrid87_90MonthlyDat/TAUXTAUYSSH/CESM_regrid_TAUXTAUYSSH_87_01-90_12.nc'; 
% infile1 = '../raw_data/CESMncFilesGlobalFromAgartha_ZhangQiuying/CESM_regrid_MonthlyOutput_TAUXTAUYSSH_88-01.nc'; 
infile2 = '../raw_data/CESM_GlobalRegrid87_90MonthlyDat/TAUXTAUYSSH/CESM_regrid_WVEL_87_01-90_12.nc'; 

% ncdisp(infile1)
outfile1 = ['../data_after_manipulation/LINLEkmanWvelGlobal_',date_str,'.mat'];

pic1 = ['../pics/TauAmpCurlRoContoursGlobalBorealSummerWinter_',date_str,'.png'];
pic2 = ['../pics/WvelContours_GlobalBorealSummerWinter_',date_str,'.png'];

  addpath(genpath('Func4EkmanProject/'))
% addpath(genpath('/home/dapengli@ad.geos.tamu.edu/MatlabCodes4DatAnalysis_DL'))

Constants4CESM
%=============================


%% === load data ===
disp('loading data')
% vardata = ncread(source,varname,start,count,stride)
% lon_1d_raw = ncread(infile1,'lon');
% lat_1d_raw = ncread(infile1,'lat');
taux_raw = ncread(infile1,'TAUX',start_3dvar,count_3dvar,stride_3dvar).*0.1; 
tauy_raw = ncread(infile1,'TAUY',start_3dvar,count_3dvar,stride_3dvar).*0.1; 
ssh_raw  = ncread(infile1,'SSH', start_3dvar,count_3dvar,stride_3dvar).*0.01; % 
% wvel_raw = ncread(infile2,'WVEL',start_4dvar,count_4dvar,stride_4dvar).*0.01; 
% ==================


%% === data analysis ===
% --- compute geostrophic velocity and geostrophic vorticity ---
disp('computing geostrophic vel and vorticity')
for iday = 1 : length(jultime_raw) % loop through time
  % compute geostrophic velocity
  [Ug_raw(:,:,iday),Vg_raw(:,:,iday)] = CalcGeostrophyVel4UnevenGridsFunc( ...
       ssh_raw(:,:,iday)',x_2d,y_2d,lat_1d); 
    
  % compute geostrophic vorticity
  kesai_raw(:,:,iday) = CalcCurlz4UnevenGridsFunc(x_2d,y_2d, ...
           Ug_raw(:,:,iday),Vg_raw(:,:,iday));
       
  % compute Rossby num
    Ro_raw(:,:,iday) = kesai_raw(:,:,iday)./f_rp; 
end
% --------------------------------------------

% --- compute wind stress amplitude and curl ---
disp('compute wind stress amplitude and curl')
for iday = 1 : length(jultime_raw) % loop through time
  % compute wind stress amplitude 
  tau_amp_raw(:,:,iday) = sqrt((taux_raw(:,:,iday)').^2+(tauy_raw(:,:,iday)').^2);   
  % compute wind stress curl
  tau_curl_raw(:,:,iday)= CalcCurlz4UnevenGridsFunc(x_2d,y_2d, ...
     taux_raw(:,:,iday)',tauy_raw(:,:,iday)');
end
% ----------------------------------------------

% --- compute linear and nonlinear Ekman Wvel ---
disp('compute Ekman Wvel')
for iday = 1 : length(jultime_raw) % loop through time
% compute linear Ekman Wvel
  [w_nl_raw(:,:,iday)] = CalcEkmanWvelFunc(rho_w,x_2d,y_2d, ...
     taux_raw(:,:,iday)',tauy_raw(:,:,iday)',f_rp,kesai_raw(:,:,iday),f_min);
% compute nonlinear Ekman Wvel 
  [w_li_raw(:,:,iday)] = CalcEkmanWvelFunc(rho_w,x_2d,y_2d, ...
     taux_raw(:,:,iday)',tauy_raw(:,:,iday)',f_rp,0,f_min);
end

w_df_raw = w_nl_raw - w_li_raw; % df: difference between high and low CESM model results
% -------------------------

% --- compute summer and winter mean ---
% --- winter ---  
jultime_wt = jultime_raw(IndxSea{1});
tau_amp_wt = tau_amp_raw(:,:,IndxSea{1});
tau_curl_wt = tau_curl_raw(:,:,IndxSea{1});
% kesai_wt = kesai_raw(:,:,IndxSea{1});
Ro_wt  = Ro_raw(:,:,IndxSea{1});
w_nl_wt = w_nl_raw(:,:,IndxSea{1});
w_li_wt = w_li_raw(:,:,IndxSea{1});
w_df_wt = w_df_raw(:,:,IndxSea{1});
% w_ce_wt = w_ce_raw(:,:,IndxSea{1});
% ---------------

% --- summer ---  
jultime_su = jultime(IndxSea{2});
tau_amp_su = tau_amp_raw(:,:,IndxSea{2});
tau_curl_su = tau_curl_raw(:,:,IndxSea{2});
% kesai_su = kesai_raw(:,:,IndxSea{2});
Ro_su  = Ro_raw(:,:,IndxSea{2});
w_nl_su = w_nl_raw(:,:,IndxSea{2});
w_li_su = w_li_raw(:,:,IndxSea{2});
% w_ce_wt = w_ce_raw(:,:,IndxSea{2});
w_df_su = w_df_raw(:,:,IndxSea{2});
% ---------------
% --------------------------------------
% ======================

  
%% === make pics ===
f1=figure; % wind stress amp and curl

   set(f1,'units','normalized','position',[0,0,1,1])
   
subplot(2,2,1)
  pcolor(lon_1d,lat_1d,tau_amp_su);shading interp;
   
subplot(2,2,2)
  pcolor(lon_1d,lat_1d,tau_amp_wt);shading interp;
  
subplot(2,2,3)
  pcolor(lon_1d,lat_1d,tau_curl_su);shading interp;
   
subplot(2,2,4)
  pcolor(lon_1d,lat_1d,tau_curl_wt);shading interp;
  
  
% figure;
%    cbar_limits_Wvel = [-5 5].*1e-5;
% %  set(gca,'color','k');
%  ax=subplot(1,1,1);%set(ax,'color','k');
% %  cbar_limits = [-0.5 0.5] % for wind stress
%  cbar_limits_Ro = [-0.4 0.4]; % for Ro 
%  
%  hold on;
% %  m_proj('miller','lon',[0 360],'lat',[-70 70]);
%  m_proj('miller','lon',lon_limits,'lat',lat_limits);
%  %m_grid('linestyle','-','gridcolor','w');
% 
% %  m_pcolor(lon_1d,lat_1d,Ro_raw(:,:,1));shading interp; 
%    m_pcolor(lon_1d,lat_1d,w_nl_raw(:,:,1)); shading interp;
%  caxis([cbar_limits_Wvel]);colorbar;
%  m_coast('patch',[.7 .7 .7],'edgecolor','none');
%   polarmap;
%   
%    m_grid('linestyle','none')
%===================


% %% === make pics ===
% f2=figure; 
% 
%   pic2_size=[10 10]; 
% % pic size: unit [inch] 3.155 inch is the width(length in x direction) of JPO 1 column pic, height varies according to your needs
%   font_size=12;
%   set(f2,'Units','inches','Position',[5,5,pic2_size]);
% % 'Position',[5,5,pic1_size] unit is inch, 5,5, is the left-bottom corner
% % of the fig on screen(no need to chage), pic_size specify the size of the window 
% % as the same size of the output fig.
%   W_lim = [-2e-5 2e-5];
% 
% % ~~~ generate subplot position ~~~  
%   row_num=3;col_num=2;margin_left=0.06;
%   margin_right=0.03;margin_top=0.04;margin_botm=0.06;
%   pics_dist_x=0.03; pics_dist_y=0.05;
%  
%   [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
%       margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
%  
% ax1=subplot('Position',sbplt_posit(1,:));
%   pcolor(Lon,Lat,tau_st_av);shading interp;
%   caxis([0 0.5]);colormap(ax1,polarmap(64,1));hc=colorbar;title(hc,'N/m^2');
%   set(gca,'xlim',[130 170],'xticklabel',[],'ylim',[26.5 45.5],'fontsize',font_size);
%   ylabel('Lat','fontsize',font_size);
%   title([num2str(length(indx_st)),' storms composite-av \tau']);
% 
% ax2=subplot('Position',sbplt_posit(2,:));
%   pcolor(Lon,Lat,Ro_st_av);shading interp;
%   caxis([-0.3 0.3]);colormap(ax2,polarmap(64,1));
%     hcb=colorbar;set(hcb,'ylim',[-0.3 0.3],'YTick',[-0.3:0.1:0.3])
%     hold on;
%   [c,h]=contour(Lon,Lat,Ro_st_std,[Ro_std_limit Ro_std_limit],'k');
%   set(h,'LineWidth',1.5);
%   set(gca,'xlim',[130 170],'ylim',[26.5 45.5],'xticklabel',[], ...
%       'yticklabel',[],'fontsize',font_size);
%   title(['Ro']);
% 
% ax3=subplot('Position',sbplt_posit(3,:));
%   pcolor(Lon,Lat,w_li_st_av);shading interp;
%   caxis(W_lim);colormap(ax3,polarmap(64,1));hc=colorbar;title(hc,'m/s');
% %     hold on;
% %   [c,h]=contour(Lon,Lat,Ro_st_std,[Ro_std_limit Ro_std_limit],'k');
% %   set(h,'LineWidth',1.5);
%   set(gca,'xlim',[130 170],'xticklabel',[],'ylim',[26.5 45.5],'fontsize',font_size);
%     ylabel('Lat','fontsize',font_size);title(['W_f']);  
%   
% ax4=subplot('Position',sbplt_posit(4,:));
%   pcolor(Lon,Lat,w_nl_st_av);shading interp;
% %     hold on;
% %   [c,h]=contour(Lon,Lat,Ro_st_std,[Ro_std_limit Ro_std_limit],'k');
% %   set(h,'LineWidth',1.5);
%   caxis(W_lim);colormap(ax4,polarmap(64,1));hc=colorbar;title(hc,'m/s');
%   set(gca,'xlim',[130 170],'xticklabel',[],'ylim',[26.5 45.5],'yticklabel',[],'fontsize',font_size);
%   title(['W_{\xi}']);
%   
% ax5=subplot('Position',sbplt_posit(5,:));
%   pcolor(Lon,Lat,w_df_st_av);shading interp;
%     hold on;
%   [c,h]=contour(Lon,Lat,Ro_st_std,[Ro_std_limit Ro_std_limit],'k');
%   set(h,'LineWidth',1.5);
%   caxis(W_lim);colormap(ax5,polarmap(64,1));hc=colorbar;title(hc,'m/s');
%   set(gca,'xlim',[130 170],'ylim',[26.5 45.5],'fontsize',font_size);
%     xlabel('Lon','fontsize',font_size);ylabel('Lat','fontsize',font_size);
%     title(['W_{\xi}-W_f']);
%   
% ax6=subplot('Position',sbplt_posit(6,:));
%   pcolor(Lon,Lat,w_ce_st_av);shading interp;
%     hold on;
%   [c,h]=contour(Lon,Lat,Ro_st_std,[Ro_std_limit Ro_std_limit],'k');
%    set(h,'LineWidth',1.5);
%   caxis(W_lim);colormap(ax6,polarmap(64,1));hc=colorbar;title(hc,'m/s'); 
%   set(gca,'xlim',[130 170],'ylim',[26.5 45.5],'yticklabel',[],'fontsize',font_size);
%   xlabel('Lon','fontsize',font_size);
%   title(['W_{cesm} depth:',num2str(zW_raw(indxZ)),'m']); 
% ====================

  
%% === output data ===
print(f1,'-dpng',pic1) 
% print(f2,'-dpng','-r500',pic2) 
% ====================