%% ===readme===

% descrip: matlab scripts plot 

% update history:
% v1.0 DL 2019Nov14

% extra notes:
% =============


%% ====set up environments====
clear all;close all;clc;

  date_str='2019Dec10';
  
infile1 = '../data_after_manipulation/CESMtaux_KE_2019Dec01.mat';  
infile2 = '../data_after_manipulation/CESMtauy_KE_2019Dec01.mat';  
infile3 = '../data_after_manipulation/CESMuvel_KE_2019Dec01.mat';  
infile4 = '../data_after_manipulation/CESMvvel_KE_2019Dec01.mat';
infile5 = '../data_after_manipulation/CESMwvel_KE_2019Dec01.mat';
infile6 = '../data_after_manipulation/CESMtemp_KE_2019Dec01.mat';

pic1 = ['../pics/DpthProfWvelQekman_KE_',date_str,'.png'];

addpath(genpath('Func4EkmanProject/'))
addpath(genpath('/home/dapengli@ad.geos.tamu.edu/MatlabCodes4DatAnalysis_DL'))
rho_w = 1020;

% lat and lon range
Lon_vec = [131:0.1:169];
Lat_vec = [45 :-0.1: 27]';
Ro_std_limit = 0.15; 
tau_st_limt = 0.25; 
%=============================


%% === load data ===
TAUX_CE = load(infile1);
TAUY_CE = load(infile2);
UVEL_CE = load(infile3);
VVEL_CE = load(infile4);
WVEL_CE = load(infile5);
TEMP_CE = load(infile6);
%===================


%% === data analysis ===
% convert lat and lon to x and y
[Lon,Lat]=meshgrid(Lon_vec,Lat_vec);

% Lat_vec=U_CE.Lat_vec;Lon_vec=U_CE.Lon_vec;
% Lon=U_CE.Lon;Lat=U_CE.Lat;
[x_vec,y_vec,x,y] = LatLon2XYFunc(Lon_vec,Lat_vec);
f = repmat(sw_f(Lat_vec),1,length(Lon_vec));

zU_raw   = UVEL_CE.z;
zW_raw   = WVEL_CE.z; 
zT_raw   = TEMP_CE.z;
jultime_vec = UVEL_CE.jultime_vec;
jultime= datenum(jultime_vec);

% --- winter ---
  IndxMon = FindMonthlyTimeIndxFunc(jultime_vec);
  
% % Sea: Season
%   IndxSea{1} = [IndxMon{1};IndxMon{2};IndxMon{3}; ...
%       IndxMon{10};IndxMon{11};IndxMon{12}]; % Winter season (Jing et al. 2019)
  IndxSea{1} = [IndxMon{1};IndxMon{2}; ...
      IndxMon{12}]; % Winter season (Jing et al. 2019)

jultime_wt = jultime(IndxSea{1});  
taux_wt = TAUX_CE.taux(:,:,IndxSea{1});
tauy_wt = TAUY_CE.tauy(:,:,IndxSea{1});
u_ce_wt = UVEL_CE.u(:,:,:,IndxSea{1});
v_ce_wt = VVEL_CE.v(:,:,:,IndxSea{1});
w_ce_wt = WVEL_CE.w(:,:,:,IndxSea{1});
t_ce_wt = TEMP_CE.temp(:,:,:,IndxSea{1});

clear TAUX_CE TAUY_CE UVEL_CE VVEL_CE WVEL_CE TEMP_CE

% iz=2; itime=3; a=v_ce_st(:,:,iz,itime); whos a

for itime = 1 : length(jultime_wt) % loop through time
    
    itime
    
    for iz = 1 : length(zT_raw) % loop through depth
        
      % compute relative vorticity
      kesai = CalcRelVorticityFunc(x,y, ...
         u_ce_wt(:,:,iz,itime),v_ce_wt(:,:,iz,itime));
   
      Ro_wt(:,:,iz,itime) = kesai./f; % Rossby num


     % compute linear and nonlinear Ekman Wvel
     [w_nl_wt(:,:,iz,itime)] = CalcNonLinearEkmanWvelFunc(rho_w,x,y, ...
       taux_wt(:,:,itime),tauy_wt(:,:,itime),f,kesai);
 
     [w_li_wt(:,:,iz,itime)] = CalcNonLinearEkmanWvelFunc(rho_w,x,y, ...
       taux_wt(:,:,itime),tauy_wt(:,:,itime),f,0);
   
    end
    
end

for iz = 1 : length(zT_raw) % loop through depth
    iz
  [Q_li_wt(:,:,iz,:),~,~] = CalcVerticalEddyHeatFluxFunc( ...
      squeeze(w_li_wt(:,:,iz,:)),squeeze(t_ce_wt(:,:,iz,:)),rho_w);
  % Q_nl_tmav_wt = nanmean(Q_nl_wt,3);
  % Q_nl_llav_wt = squeeze(nanmean(nanmean(Q_nl_wt,1),2));


  [Q_nl_wt(:,:,iz,:),~,~] = CalcVerticalEddyHeatFluxFunc( ...
      squeeze(w_nl_wt(:,:,iz,:)),squeeze(t_ce_wt(:,:,iz,:)),rho_w);
  % Q_nl_tmav_wt = nanmean(Q_nl_wt,3);
  % Q_nl_llav_wt = squeeze(nanmean(nanmean(Q_nl_wt,1),2));

  [Q_ce_wt(:,:,iz,:),~,~]= CalcVerticalEddyHeatFluxFunc( ...
      squeeze(w_ce_wt(:,:,iz,:)),squeeze(t_ce_wt(:,:,iz,:)),rho_w);
  % Q_ce_tmav_wt = nanmean(Q_ce_wt,3);
  % Q_ce_llav_wt = squeeze(nanmean(nanmean(Q_ce_wt,1),2));

end

% --- composite av during winter storm events ---
tau_wt  = sqrt(taux_wt.^2+tauy_wt.^2);
tau_llav_wt = squeeze(nanmean(nanmean(tau_wt,1),2)); % llav: lat and lon av
indx_st = find(tau_llav_wt>=tau_st_limt); % length(indx_st) = 76

% taux_st = taux_wt(:,:,indx_st);
% tauy_st = tauy_wt(:,:,indx_st);
% u_ce_st = u_ce_wt(:,:,:,indx_st);
% v_ce_st = v_ce_wt(:,:,:,indx_st);
% w_ce_st = w_ce_wt(:,:,:,indx_st);
% t_ce_st = t_ce_wt(:,:,:,indx_st);

clear taux_wt tauy_wt u_ce_wt v_ce_wt t_ce_wt

Ro_st = Ro_wt(:,:,:,indx_st);
Ro_st_timeav = nanmean(Ro_st,4);
% Ro_st_std = nanstd(Ro_st,0,3);

w_li_st = w_li_wt(:,:,:,indx_st);
w_li_st_timeav = nanmean(w_li_st,4);
w_nl_st = w_nl_wt(:,:,:,indx_st);
w_nl_st_timeav = nanmean(w_nl_st,4);

% w_df_st = w_df_wt(:,:,indx_st);
% w_df_st_av = nanmean(w_df_st,3);
w_ce_st = w_ce_wt(:,:,:,indx_st);
w_ce_st_timeav = nanmean(w_ce_st,4);

% kesai_st = kesai_wt(:,:,indx_st);
% kesai_std_st = nanstd(kesai_st,0,3);

Q_li_st = Q_li_wt(:,:,:,indx_st);
Q_li_st_timeav = nanmean(Q_li_st,4);

Q_nl_st = Q_nl_wt(:,:,:,indx_st);
Q_nl_st_timeav = nanmean(Q_nl_st,4);

% Q_df_st = Q_nl_st - Q_li_st;
% Q_df_st_av = nanmean(Q_df_st,3);

Q_ce_st = Q_ce_wt(:,:,:,indx_st);
Q_ce_st_timeav = nanmean(Q_ce_st,4);
% -------------------

indxZ = 3;
zU_raw(indxZ),zW_raw(indxZ),

Ro_st_std = nanstd(squeeze(Ro_st(:,:,indxZ,:)),0,3);

% figure % compare it with WvelContours_KE_2019Dec05.png Ro subplot
%   pcolor(Lon,Lat,Ro_st_timeav(:,:,indxZ));shading interp;
%   caxis([-0.3 0.3]);polarmap(64,1);colorbar;hold on;
%   [c,h]=contour(Lon,Lat,Ro_st_std,[Ro_std_limit Ro_std_limit],'k');


% --- Ro > Ro std ---
indxRo = find(Ro_st_std>Ro_std_limit); % length(indxRo) = 6884;

for iz = 1 : length(zT_raw)
    
    iz
    
    Ro_st_timeav_iz = Ro_st_timeav(:,:,iz);
    [Ro_st_z_min(iz),Ro_st_z_av(iz),Ro_st_z_max(iz)]=bootstrap5( ...
        Ro_st_timeav_iz(indxRo));
    
    w_li_st_timeav_iz = w_li_st_timeav(:,:,iz);
    [w_li_st_z_min(iz),w_li_st_z_av(iz),w_li_st_z_max(iz)]=bootstrap5( ...
        w_li_st_timeav_iz(indxRo));
    
    w_nl_st_timeav_iz = w_nl_st_timeav(:,:,iz);
    [w_nl_st_z_min(iz),w_nl_st_z_av(iz),w_nl_st_z_max(iz)]=bootstrap5( ...
        w_nl_st_timeav_iz(indxRo));
    
    w_ce_st_timeav_iz = w_ce_st_timeav(:,:,iz);
    [w_ce_st_z_min(iz),w_ce_st_z_av(iz),w_ce_st_z_max(iz)]=bootstrap5( ...
        w_ce_st_timeav_iz(indxRo));
    
    Q_li_st_timeav_iz = Q_li_st_timeav(:,:,iz);
    [Q_li_st_z_min(iz),Q_li_st_z_av(iz),Q_li_st_z_max(iz)]=bootstrap5( ...
        Q_li_st_timeav_iz(indxRo));
    
    Q_nl_st_timeav_iz = Q_nl_st_timeav(:,:,iz);
    [Q_nl_st_z_min(iz),Q_nl_st_z_av(iz),Q_nl_st_z_max(iz)]=bootstrap5( ...
        Q_nl_st_timeav_iz(indxRo));
    
    Q_ce_st_timeav_iz = Q_ce_st_timeav(:,:,iz);
    [Q_ce_st_z_min(iz),Q_ce_st_z_av(iz),Q_ce_st_z_max(iz)]=bootstrap5( ...
        Q_ce_st_timeav_iz(indxRo));
    
    clear Ro_st_timeav_iz w_li_st_timeav_iz w_nl_st_timeav_iz w_ce_st_timeav_iz
    clear Q_li_st_timeav_iz Q_nl_st_timeav_iz Q_ce_st_timeav_iz
    
end
% ===================


%% === make pics ===
f1=figure; % composite av
 
  pic1_size=[8 4]; % pic size: unit [inch] 
  font_size=12;
  set(f1,'Units','inches','Position',[5,5,pic1_size]);
% 'Position',[5,5,pic1_size] unit is inch, 5,5, is the left-bottom corner
% of the fig on screen(no need to chage), pic_size specify the size of the window 
% as the same size of the output fig. 
 
% ~~~ generate subplot position ~~~  
  row_num=1;col_num=3;margin_left=0.1;
  margin_right=0.05;margin_top=0.1;margin_botm=0.15;
  pics_dist_x=0.05; pics_dist_y=0.0;
 
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 
 subplot('Position',sbplt_posit(1,:));
  [h_fil h_plt]=plt_shade_line(Ro_st_z_min',Ro_st_z_av', ...
      Ro_st_z_max',-zU_raw,'b');grid on;
    set(h_plt,'linewidth',1.5);
    set(gca,'xlim',[0.01 0.03],'xtick',[0.01:0.01:0.03],'ylim',[-80 -20], ...
        'ytick',[-75:10:-25],'fontsize',font_size);
    xlabel('Ro','fontsize',font_size);ylabel('Depth [m]','fontsize',font_size);
    title('Rossby Num Prof');
    
 subplot('Position',sbplt_posit(2,:));
   [h_fil h_plt1]=plt_shade_line(w_li_st_z_min',w_li_st_z_av', ...
       w_li_st_z_max',-zU_raw,'b');
   set(h_plt1,'linewidth',1.5);
   [h_fil h_plt2]=plt_shade_line(w_nl_st_z_min',w_nl_st_z_av', ...
       w_nl_st_z_max',-zU_raw,'r');
   set(h_plt2,'linewidth',1.5);
   [h_fil h_plt3]=plt_shade_line(w_ce_st_z_min',w_ce_st_z_av', ...
       w_ce_st_z_max',-zU_raw,'k');
   set(h_plt3,'linewidth',1.5);
    legend([h_plt1,h_plt2,h_plt3],{'W_f','W_{\xi}','W_{cesm}'}, ...
       'fontsize',9,'Position',[0.53 0.69 0.11 0.17]);legend boxoff;
    set(gca,'xlim',[0 6e-6],'xtick',[0:2:6].*1e-6,'ylim',[-80 -20], ...
        'ytick',[-75:10:-25],'yticklabel',[],'fontsize',font_size);grid on;
    xlabel('W [m/s]','fontsize',font_size);%ylabel('Depth [m]');
    title('W Dpth Prof');
    
 subplot('Position',sbplt_posit(3,:));
   [h_fil h_plt1]=plt_shade_line(Q_li_st_z_min',Q_li_st_z_av', ...
       Q_li_st_z_max',-zU_raw,'b');grid on;
    set(h_plt1,'linewidth',1.5);
   [h_fil h_plt2]=plt_shade_line(Q_nl_st_z_min',Q_nl_st_z_av', ...
       Q_nl_st_z_max',-zU_raw,'r');
    set(h_plt2,'linewidth',1.5);
   [h_fil h_plt3]=plt_shade_line(Q_ce_st_z_min',Q_ce_st_z_av', ...
       Q_ce_st_z_max',-zU_raw,'k');
    set(h_plt3,'linewidth',1.5);
    legend([h_plt1,h_plt2,h_plt3],{'Q_f','Q_{\xi}','Q_{cesm}'}, ...
        'fontsize',9,'Position',[0.77 0.21 0.11 0.17]);legend boxoff;
    set(gca,'xlim',[0 270],'xtick',[0:50:250],'ylim',[-80 -20], ...
        'ytick',[-75:10:-25],'yticklabel',[],'fontsize',font_size);
    xlabel('Q [W/m^2]','fontsize',font_size);%ylabel('Depth [m]');
    title('Vertical Heat Flux Prof')
% ======================

  
%% === output data ===
print(f1,'-dpng','-r500',pic1)
% ====================