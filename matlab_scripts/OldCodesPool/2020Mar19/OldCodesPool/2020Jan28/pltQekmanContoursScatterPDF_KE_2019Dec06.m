%% ===readme===

% descrip: matlab scripts plot tau, UVW vel, Kesai from 
% 1) CESM model output
% 2) Stern nonlinear Ekman vertical pumping velocity
% 3) linear Ekman vertical pumping velocity

% update history:
% v1.0 DL 2019Nov14

% extra notes:
% =============


%% ====set up environments====
clear all;close all;clc;

  date_str='2019Dec05';
  
infile1 = '../data_after_manipulation/CESMtaux_KE_2019Dec01.mat';  
infile2 = '../data_after_manipulation/CESMtauy_KE_2019Dec01.mat';  
infile3 = '../data_after_manipulation/CESMuvel_KE_2019Dec01.mat';  
infile4 = '../data_after_manipulation/CESMvvel_KE_2019Dec01.mat';
infile5 = '../data_after_manipulation/CESMwvel_KE_2019Dec01.mat';
infile6 = '../data_after_manipulation/CESMtemp_KE_2019Dec01.mat';

pic1 = ['../pics/QekmanContours_KE_',date_str,'.png'];
pic2 = ['../pics/QekmanScatter_KE_',date_str,'.png'];
pic3 = ['../pics/QekmanPDF_KE_',date_str,'.png'];

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

zU_raw   = UVEL_CE.z;
zW_raw   = WVEL_CE.z; 
zT_raw   = TEMP_CE.z;
indxZ = 3;
zU_raw(indxZ),zW_raw(indxZ),zT_raw(indxZ),

u_ce_raw = squeeze(UVEL_CE.u(:,:,indxZ,:));
v_ce_raw = squeeze(VVEL_CE.v(:,:,indxZ,:));
w_ce_raw = squeeze(WVEL_CE.w(:,:,indxZ,:));
t_ce_raw = squeeze(TEMP_CE.temp(:,:,indxZ,:));

taux_raw = TAUX_CE.taux;
tauy_raw = TAUY_CE.tauy;
tau_raw  = sqrt(taux_raw.^2+tauy_raw.^2);
tau_llav = squeeze(nanmean(nanmean(tau_raw,1),2)); % llav: lat and lon av

jultime_vec = UVEL_CE.jultime_vec;
jultime= datenum(jultime_vec);

clear TAUX_CE TAUY_CE UVEL_CE VVEL_CE WVEL_CE TEMP_CE
f = repmat(sw_f(Lat_vec),1,length(Lon_vec));

% compute relative vorticity
for i = 1 : length(jultime) % loop through time
    kesai_raw(:,:,i) = CalcRelVorticityFunc(x,y, ...
       u_ce_raw(:,:,i),v_ce_raw(:,:,i));
   
    Ro_raw(:,:,i) = kesai_raw(:,:,i)./f; % Rossby num
end

% compute linear and nonlinear Ekman Wvel
for i = 1 : length(jultime) % loop through time
    
 [w_nl_raw(:,:,i)] = CalcNonLinearEkmanWvelFunc(rho_w,x,y, ...
     taux_raw(:,:,i),tauy_raw(:,:,i),f,kesai_raw(:,:,i));
 
 [w_li_raw(:,:,i)] = CalcNonLinearEkmanWvelFunc(rho_w,x,y, ...
     taux_raw(:,:,i),tauy_raw(:,:,i),f,0);
 
end

w_df_raw = w_nl_raw - w_li_raw; % df: difference between high and low CESM model

%%

% find specific time (winter season)
  IndxMon = FindMonthlyTimeIndxFunc(jultime_vec);
  
% Sea: Season
  IndxSea{1} = [IndxMon{1};IndxMon{2};IndxMon{3}; ...
      IndxMon{10};IndxMon{11};IndxMon{12}]; % Winter season (Jing et al. 2019)
  IndxSea{2} = [IndxMon{4};IndxMon{5};IndxMon{6}; ...
      IndxMon{7};IndxMon{8};IndxMon{9}]; % Winter season (Jing et al. 2019)

%   IndxSea{1} = [IndxMon{1};IndxMon{2}; ...
%       IndxMon{12}]; % Winter season (Jing et al. 2019)
  
% --- winter ---  
jultime_wt = jultime(IndxSea{1});

tau_llav_wt  = tau_llav(IndxSea{1});
tau_wt       = tau_raw(:,:,IndxSea{1});

kesai_wt = kesai_raw(:,:,IndxSea{1});
Ro_wt  = Ro_raw(:,:,IndxSea{1});

w_nl_wt = w_nl_raw(:,:,IndxSea{1});
w_li_wt = w_li_raw(:,:,IndxSea{1});
w_df_wt = w_df_raw(:,:,IndxSea{1});
w_ce_wt = w_ce_raw(:,:,IndxSea{1});
t_ce_wt = t_ce_raw(:,:,IndxSea{1});

[Q_li_wt,~,~] = CalcVerticalEddyHeatFluxFunc(w_li_wt,t_ce_wt,rho_w);
% Q_nl_tmav_wt = nanmean(Q_nl_wt,3);
% Q_nl_llav_wt = squeeze(nanmean(nanmean(Q_nl_wt,1),2));


[Q_nl_wt,~,~] = CalcVerticalEddyHeatFluxFunc(w_nl_wt,t_ce_wt,rho_w);
% Q_nl_tmav_wt = nanmean(Q_nl_wt,3);
% Q_nl_llav_wt = squeeze(nanmean(nanmean(Q_nl_wt,1),2));

[Q_ce_wt,~,~]= CalcVerticalEddyHeatFluxFunc(w_ce_wt,t_ce_wt,rho_w);
% Q_ce_tmav_wt = nanmean(Q_ce_wt,3);
% Q_ce_llav_wt = squeeze(nanmean(nanmean(Q_ce_wt,1),2));
% ---------------

% --- summer ---  
jultime_sm = jultime(IndxSea{2});

tau_llav_sm  = tau_llav(IndxSea{2});
tau_sm       = tau_raw(:,:,IndxSea{2});

kesai_sm = kesai_raw(:,:,IndxSea{2});
Ro_sm  = Ro_raw(:,:,IndxSea{2});

w_nl_sm = w_nl_raw(:,:,IndxSea{2});
w_li_sm = w_li_raw(:,:,IndxSea{2});
w_df_sm = w_df_raw(:,:,IndxSea{2});
w_ce_sm = w_ce_raw(:,:,IndxSea{2});
t_ce_sm = t_ce_raw(:,:,IndxSea{2});

[Q_li_sm,~,~] = CalcVerticalEddyHeatFluxFunc(w_li_sm,t_ce_sm,rho_w);
% Q_nl_tmav_sm = nanmean(Q_nl_sm,3);
% Q_nl_llav_sm = squeeze(nanmean(nanmean(Q_nl_sm,1),2));


[Q_nl_sm,~,~] = CalcVerticalEddyHeatFluxFunc(w_nl_sm,t_ce_sm,rho_w);
% Q_nl_tmav_sm = nanmean(Q_nl_sm,3);
% Q_nl_llav_sm = squeeze(nanmean(nanmean(Q_nl_sm,1),2));

[Q_ce_sm,~,~]= CalcVerticalEddyHeatFluxFunc(w_ce_sm,t_ce_sm,rho_w);
% Q_ce_tmav_sm = nanmean(Q_ce_sm,3);
% Q_ce_llav_sm = squeeze(nanmean(nanmean(Q_ce_sm,1),2));
% --------------------


% ---------------

% composite av during winter storm events
indx_st = find(tau_llav_wt>=tau_st_limt); % length(indx_st) = 6884

tau_st = tau_wt(:,:,indx_st);
tau_st_av = nanmean(tau_st,3);

Ro_st = Ro_wt(:,:,indx_st);
Ro_st_av = nanmean(Ro_st,3);
Ro_st_std = nanstd(Ro_st,0,3);

w_li_st = w_li_wt(:,:,indx_st);
w_li_st_av = nanmean(w_li_st,3);
w_nl_st = w_nl_wt(:,:,indx_st);
w_nl_st_av = nanmean(w_nl_st,3);

% w_df_st = w_df_wt(:,:,indx_st);
% w_df_st_av = nanmean(w_df_st,3);
w_ce_st = w_ce_wt(:,:,indx_st);
w_ce_st_av = nanmean(w_ce_st,3);

% kesai_st = kesai_wt(:,:,indx_st);
% kesai_std_st = nanstd(kesai_st,0,3);

Q_li_st = Q_li_wt(:,:,indx_st);
Q_li_st_av = nanmean(Q_li_st,3);

Q_nl_st = Q_nl_wt(:,:,indx_st);
Q_nl_st_av = nanmean(Q_nl_st,3);

Q_df_st = Q_nl_st - Q_li_st;
Q_df_st_av = nanmean(Q_df_st,3);

Q_ce_st = Q_ce_wt(:,:,indx_st);
Q_ce_st_av = nanmean(Q_ce_st,3);

% --- Ro > Ro std ---
indxRo = find(Ro_st_std>Ro_std_limit);length(indxRo)
w_li_Ro = w_li_st_av(indxRo);
w_nl_Ro = w_nl_st_av(indxRo); 
w_ce_Ro = w_ce_st_av(indxRo); 
Q_li_Ro = Q_li_st_av(indxRo);
Q_nl_Ro = Q_nl_st_av(indxRo); 
Q_ce_Ro = Q_ce_st_av(indxRo); 

% xi=w_nl_Ro;yi=w_li_Ro;x0=0;y0=0;m=1;
% % P = polyfix(xi,yi,x0,y0,m)
% P_nl2li = polyfix(xi(~isnan(yi)),yi(~isnan(yi)),x0,y0,m); 

% brob = robustfit(X,Y), brob(1)+brob(2)*x
bW_nl2li = robustfit(w_nl_Ro,w_li_Ro); % y = -2.6938e-07 + 0.1786*x
bW_nl2ce = robustfit(w_nl_Ro,w_ce_Ro); % y = 8.0752e-07  + 0.4986*x
bQ_nl2li = robustfit(Q_nl_Ro,Q_li_Ro); % y = 17.0900    + 0.1648*x
bQ_nl2ce = robustfit(Q_nl_Ro,Q_ce_Ro); % y = 160.6539   + 0.9618*x

indxNoNaN_Wnl = find(~isnan(w_nl_Ro));
indxNoNaN_Wli = find(~isnan(w_li_Ro));
indxNoNaN_Wce = find(~isnan(w_ce_Ro));
% union([1 2 3], [2 3 4 5 6])
indxNoNaN_W = intersect(intersect(indxNoNaN_Wnl,indxNoNaN_Wli),indxNoNaN_Wce);

indxNoNaN_Qnl = find(~isnan(Q_nl_Ro));
indxNoNaN_Qli = find(~isnan(Q_li_Ro));
indxNoNaN_Qce = find(~isnan(Q_ce_Ro));
indxNoNaN_Q = intersect(intersect(indxNoNaN_Qnl,indxNoNaN_Qli),indxNoNaN_Qce);

rW_nl2li = corr2(w_nl_Ro(indxNoNaN_W),w_li_Ro(indxNoNaN_W)); % rW_nl2li=0.5174
rW_nl2ce = corr2(w_nl_Ro(indxNoNaN_W),w_ce_Ro(indxNoNaN_W)); % rW_nl2ce=0.1728
rQ_nl2li = corr2(Q_nl_Ro(indxNoNaN_Q),Q_li_Ro(indxNoNaN_Q)); % rQ_nl2li=0.5223
rQ_nl2ce = corr2(Q_nl_Ro(indxNoNaN_Q),Q_ce_Ro(indxNoNaN_Q)); % rQ_nl2ce=0.3341

% nanmean(Q_li_st_av(indxRo)) = 23.3479
% nanmean(Q_nl_st_av(indxRo)) = 40.9070
% nanmean(Q_ce_st_av(indxRo)) = 212.5425

% 40.9070./23.3479 = 1.7521
% 40.9070./212.5425 = 0.1925  

% R = corr2(Q_nl_st_av(Ro_st_std>Ro_std_limit),Q_ce_st_av(Ro_st_std>Ro_std_limit));

% PDF
WbinsCenter = [-1:0.5:1].*1e-5; 

[N_w_li_Ro,C_w_li_Ro] = hist(w_li_Ro(:),WbinsCenter); 
P_w_li_Ro = N_w_li_Ro./sum(N_w_li_Ro);

[N_w_nl_Ro,C_w_nl_Ro] = hist(w_nl_Ro(:),WbinsCenter); 
P_w_nl_Ro = N_w_nl_Ro./sum(N_w_nl_Ro);

[N_w_ce_Ro,C_w_ce_Ro] = hist(w_ce_Ro(:),WbinsCenter); %[-25:50:600]
P_w_ce_Ro = N_w_ce_Ro./sum(N_w_ce_Ro);

QbinsCenter = [-50:50:200]; 

[N_Q_li_Ro,C_Q_li_Ro] = hist(Q_li_Ro(:),QbinsCenter); 
P_Q_li_Ro = N_Q_li_Ro./sum(N_Q_li_Ro);

[N_Q_nl_Ro,C_Q_nl_Ro] = hist(Q_nl_Ro(:),QbinsCenter); 
P_Q_nl_Ro = N_Q_nl_Ro./sum(N_Q_nl_Ro);

[N_Q_ce_Ro,C_Q_ce_Ro] = hist(Q_ce_Ro(:),QbinsCenter); %[-25:50:600]
P_Q_ce_Ro = N_Q_ce_Ro./sum(N_Q_ce_Ro);
% ==================


%% === make pics ===
f1=figure('Renderer','painters'); % composite av
 
  % set pic size
  pic1_size=[10 10]; % pic size: unit [inch] 3.155 inch is the width(length in x direction) of JPO 1 column pic, height varies according to your needs
  set(f1,'Units','inches','Position',[5,5,pic1_size]);
  font_size=12;
  Q_lim = [-200 200];
 
% ~~~ generate subplot position ~~~  
  row_num=2;col_num=2;margin_left=0.06;
  margin_right=0.03;margin_top=0.04;margin_botm=0.06;
  pics_dist_x=0.03; pics_dist_y=0.06;
 
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
 subplot('position',sbplt_posit(1,:));
  pcolor(Lon,Lat,Q_li_st_av);shading interp;
    caxis(Q_lim);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
  hold on;
  [c,h]=contour(Lon,Lat,Ro_st_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',[130 170],'xticklabel',[],'ylim',[26.5 45.5], ...
        'ytick',[30:5:45],'fontsize',font_size);
    ylabel('Lat','fontsize',font_size);
    title(['Q_f  depth: ',num2str(zT_raw(indxZ)),' m'])
    
 subplot('position',sbplt_posit(2,:));
  pcolor(Lon,Lat,Q_nl_st_av);shading interp;
    caxis(Q_lim);polarmap(64,1);hc=colorbar;title(hc,'W/m^2'); 
      hold on;
  [c,h]=contour(Lon,Lat,Ro_st_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',[130 170],'xticklabel',[],'ylim',[26.5 45.5], ...
        'ytick',[30:5:45],'yticklabel',[],'fontsize',font_size);
    title(['Q_{\xi}  depth: ',num2str(zT_raw(indxZ)),' m'])
    
 subplot('position',sbplt_posit(3,:))
  pcolor(Lon,Lat,Q_df_st_av);shading interp;
    caxis(Q_lim);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');  
      hold on;
  [c,h]=contour(Lon,Lat,Ro_st_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',[130 170],'ylim',[26.5 45.5],'ytick',[30:5:45],'fontsize',font_size);
    xlabel('Lon','fontsize',font_size);
    ylabel('Lat','fontsize',font_size);
    title(['Q_{\xi} - Q_f  depth: ',num2str(zT_raw(indxZ)),' m'])
  
 subplot('position',sbplt_posit(4,:))
  pcolor(Lon,Lat,Q_ce_st_av);shading interp;
    caxis(Q_lim.*2);polarmap(64,1);hc=colorbar;title(hc,'W/m^2'); 
      hold on;
  [c,h]=contour(Lon,Lat,Ro_st_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',[130 170],'ylim',[26.5 45.5],'ytick',[30:5:45], ...
        'yticklabel',[],'fontsize',font_size);
    xlabel('Lon','fontsize',font_size);
    title(['Q_{cesm}  depth: ',num2str(zT_raw(indxZ)),' m'])

%   for itime = 1 : 5
%     indxCor = find(Q_nl_st_av>(itime-1)*50);
%     numindxCor(itime) = length(indxCor);
%     R(itime) = corr2(Q_nl_st_av(indxCor),Q_ce_st_av(indxCor));
%     clear indxCor
%   end
%   
%   figure;
%    plot([100:100:500],R)  

%%
f2 = figure('Renderer','painters'); % scatter plot
 
   % set pic size
  pic2_size=[10 7]; % pic size: unit [inch] 3.155 inch is the width(length in x direction) of JPO 1 column pic, height varies according to your needs
  set(f2,'Units','inches','Position',[5,5,pic2_size]);
  font_size=10;
  W_lim = [-0.4e-4 0.8e-4];
  Q_lim = [-200 400];
 
% ~~~ generate subplot position ~~~  
  row_num=2;col_num=2;margin_left=0.08;
  margin_right=0.03;margin_top=0.06;margin_botm=0.08;
  pics_dist_x=0.1; pics_dist_y=0.12;
 
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  
  subplot('position',sbplt_posit(1,:))
   plot(w_nl_Ro,w_li_Ro,'b*');hold on;
   plot(W_lim,bW_nl2li(1)+bW_nl2li(2).*W_lim,'r','linewidth',1.5);grid on;
    set(gca,'xlim',[-4e-5 8e-5],'ylim',[-4e-5 6e-5],'fontsize',font_size);
    xlabel('W_{\xi} [m/s]','fontsize',font_size);
    ylabel('W_f [m/s]','fontsize',font_size);
    title(['y=',num2str(bW_nl2li(1),'%.2e'),'+',num2str(bW_nl2li(2),'%.2f'), ...
        '*x, r = ', num2str(rW_nl2li,'%.2f')])
%   [x_min,x_mn,x_max]=bootstrap5(abs(w_nl_Ro./w_li_Ro))
  
  subplot('position',sbplt_posit(2,:))
   plot(w_nl_Ro,w_ce_Ro,'b*');hold on;
   plot(W_lim,bW_nl2ce(1)+bW_nl2ce(2).*W_lim,'r','linewidth',1.5);grid on;
    set(gca,'xlim',[-4e-5 8e-5],'ylim',[-4e-4 4e-4],'ytick',[-4:2:4].*1e-4, ...
        'fontsize',font_size);
    xlabel('W_{\xi} [m/s]','fontsize',font_size);
    ylabel('W_{cesm} [m/s]','fontsize',font_size);
    title(['y=',num2str(bW_nl2ce(1),'%.2e'),'+',num2str(bW_nl2ce(2),'%.2f'), ...
        '*x, r = ',num2str(rW_nl2ce,'%.2f')])

  subplot('position',sbplt_posit(3,:))
   plot(Q_nl_Ro,Q_li_Ro,'b*');hold on;
   plot(Q_lim,bQ_nl2li(1)+bQ_nl2li(2).*Q_lim,'r','linewidth',1.5);grid on;
    set(gca,'xlim',[-200 400],'ylim',[-150 150],'fontsize',font_size)
    xlabel('Q_{\xi} [W/m^2]','fontsize',font_size);
    ylabel('Q_f [W/m^2]','fontsize',font_size);
    title(['y=',num2str(bQ_nl2li(1),'%.2f'),'+',num2str(bQ_nl2li(2),'%.2f'), ...
        '*x, r = ',num2str(rQ_nl2li,'%.2f')])
    
  subplot('position',sbplt_posit(4,:))
   plot(Q_nl_Ro,Q_ce_Ro,'b*');hold on;
   plot(Q_lim,bQ_nl2ce(1)+bQ_nl2ce(2).*Q_lim,'r','linewidth',1.5);grid on;
    set(gca,'xlim',[-200 400],'ylim',[-500 1500],'fontsize',font_size)
    xlabel('Q_{\xi} [W/m^2]','fontsize',font_size);
    ylabel('Q_{cesm} [W/m^2]','fontsize',font_size);
    title(['y=',num2str(bQ_nl2ce(1),'%.2f'),'+',num2str(bQ_nl2ce(2),'%.2f'), ...
        '*x, r = ',num2str(rQ_nl2ce,'%.2f')])
    
%%    
f3 = figure('Renderer','painters'); % PDF 
  
   % set pic size
  pic3_size=[10 5]; % pic size: unit [inch] 3.155 inch is the width(length in x direction) of JPO 1 column pic, height varies according to your needs
  set(f3,'Units','inches','Position',[5,5,pic3_size]);
  font_size=10;
 
% ~~~ generate subplot position ~~~  
  row_num=2;col_num=1;margin_left=0.06;
  margin_right=0.03;margin_top=0.08;margin_botm=0.12;
  pics_dist_x=0.0; pics_dist_y=0.18;
 
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

  subplot('position',sbplt_posit(1,:))   
   hb = bar(C_w_li_Ro,[P_w_li_Ro;P_w_nl_Ro;P_w_ce_Ro]');grid on;
    set(hb(1),'FaceColor','b');set(hb(2),'FaceColor','r');
    set(hb(3),'FaceColor','k');
    legend(hb,{'W_f','W_{\xi}','W_{cesm}'},'location','northwest' ...
        );legend boxoff;
    set(gca,'fontsize',font_size);
    xlabel('W [m/s]','fontsize',font_size);
    ylabel('[%]','fontsize',font_size);
    title(['relative distribution of W, depth: ',num2str(zT_raw(indxZ)),' m']);
   
  subplot('position',sbplt_posit(2,:))
   hb = bar(C_Q_li_Ro,[P_Q_li_Ro;P_Q_nl_Ro;P_Q_ce_Ro]');grid on;
    set(hb(1),'FaceColor','b');set(hb(2),'FaceColor','r');
    set(hb(3),'FaceColor','k');
    legend(hb,{'Q_f','Q_{\xi}','Q_{cesm}'},'location','northwest');legend boxoff;
    set(gca,'fontsize',font_size);
    xlabel('Q [W/m^2]','fontsize',font_size);
    ylabel('[%]','fontsize',font_size);
    title(['relative distribution of Q, depth: ',num2str(zT_raw(indxZ)),' m']);
% ======================

  
%% === output data ===
print(f1,'-dpng','-r500',pic1)
print(f2,'-dpng','-r500',pic2) 
print(f3,'-dpng','-r500',pic3) 
% ====================