%% ===readme===

% descrip: matlab scripts plot tau, UVW vel, Kesai from 
% 1) CESM model output
% 2) Stern nonlinear Ekman vertical pumping velocity
% 3) linear Ekman vertical pumping velocity

% update history:
% v1.0 DL 2019Nov14

% extra notes:
% W: vertical velocity
% Q: vertical heat flux
% no scatter plots are made due to low correlation (~0.2) between 
% Q_nl (Nonlinear Ekman pumping heat flux) and Q_ce (total eddy heat flux)
% =============


%% ====set up environments====
clear all;close all;clc;

  date_str='2020Jan15';
  
infile1 = '../data_after_manipulation/CESMtaux_KE_2019Dec01.mat';  
infile2 = '../data_after_manipulation/CESMtauy_KE_2019Dec01.mat';  
% infile3 = '../data_after_manipulation/CESMuvel_KE_2019Dec01.mat';  
% infile4 = '../data_after_manipulation/CESMvvel_KE_2019Dec01.mat';
infile5 = '../data_after_manipulation/CESMwvel_KE_2019Dec01.mat';
infile6 = '../data_after_manipulation/CESMtemp_KE_2019Dec01.mat';
infile7 = '../data_after_manipulation/CESMssh_KE_2020Jan15.mat';

pic2 = ['../pics/QScatterWinterStorm_KE_',date_str,'.png'];
pic3 = ['../pics/QPDFWinterStorm_KE_',date_str,'.png'];

addpath(genpath('Func4EkmanProject/'))
addpath(genpath('/home/dapengli@ad.geos.tamu.edu/MatlabCodes4DatAnalysis_DL'))
rho_w = 1020;

% lat and lon range
Lon_vec = [131:0.1:169];
Lat_vec = [45 :-0.1: 27]';
Ro_std_limit = 0.15; 
% tau_st_limt = 0.20; 
% tau_ns_limt = 0.10; 
indxZ = 1;
SecPerDay = 86400;
%=============================


%% === load data ===
TAUX_CE = load(infile1);
TAUY_CE = load(infile2);
% UVEL_CE = load(infile3);
% VVEL_CE = load(infile4);
WVEL_CE = load(infile5);
TEMP_CE = load(infile6);
SSH_CE = load(infile7);
%===================


%% === data analysis ===
% --- convert lat and lon to x and y ---
[Lon,Lat]=meshgrid(Lon_vec,Lat_vec);

% Lat_vec=U_CE.Lat_vec;Lon_vec=U_CE.Lon_vec;
% Lon=U_CE.Lon;Lat=U_CE.Lat;
[x_vec,y_vec,x,y] = LatLon2XYFunc(Lon_vec,Lat_vec);
% --------------------------------------

% --- compute f and g ---
f = repmat(sw_f(Lat_vec),1,length(Lon_vec));
g = sw_g(Lat_vec(1),0); % g = sw_g(lat,z)
% -----------------------

% --- rename loaded data ---
% u_ce_raw = squeeze(UVEL_CE.u(:,:,indxZ,:));
% v_ce_raw = squeeze(VVEL_CE.v(:,:,indxZ,:));
% zU_raw   = UVEL_CE.z; zU_raw(indxZ),

jultime_vec = TAUX_CE.jultime_vec;
jultime = datenum(jultime_vec);

taux_raw = TAUX_CE.taux;
tauy_raw = TAUY_CE.tauy;
tau_raw  = sqrt(taux_raw.^2+tauy_raw.^2);
tau_llav = squeeze(nanmean(nanmean(tau_raw,1),2)); % llav: lat and lon av
eta_raw = SSH_CE.ssh;

zW_raw = WVEL_CE.z; disp(['W depth: ',num2str(zW_raw(indxZ)),'m'])
zT_raw = TEMP_CE.z; disp(['T depth: ',num2str(zT_raw(indxZ)),'m'])

w_ce_raw = squeeze(WVEL_CE.w(:,:,indxZ,:));
t_ce_raw = squeeze(TEMP_CE.temp(:,:,indxZ,:));

clear TAUX_CE TAUY_CE WVEL_CE TEMP_CE SSH_CE
% ---------------------------------

% --- compute geostrophic velocity and geostrophic vorticity ---
for i = 1 : length(jultime) % loop through time
    % compute geostrophic velocity
    [u_gs_raw(:,:,i),v_gs_raw(:,:,i)] = CalcGeostrophicVelFunc( ...
       x,y,eta_raw(:,:,i),f,g); 
    
    % compute geostrophic vorticity
    kesai_raw(:,:,i) = CalcRelVorticityFunc(x,y, ...
       u_gs_raw(:,:,i),v_gs_raw(:,:,i));
   
    Ro_raw(:,:,i) = kesai_raw(:,:,i)./f; % Rossby num
end
% --------------------------------------------

% --- compute linear and nonlinear Ekman Wvel ---
for i = 1 : length(jultime) % loop through time
    
 [w_nl_raw(:,:,i)] = CalcNonLinearEkmanWvelFunc(rho_w,x,y, ...
     taux_raw(:,:,i),tauy_raw(:,:,i),f,kesai_raw(:,:,i));
 
 [w_li_raw(:,:,i)] = CalcNonLinearEkmanWvelFunc(rho_w,x,y, ...
     taux_raw(:,:,i),tauy_raw(:,:,i),f,0);
 
end

% w_df_raw = w_nl_raw - w_li_raw; % df: difference between high and low CESM model
% -------------------------

% --- find winter season --- 
IndxMon = FindMonthlyTimeIndxFunc(jultime_vec);
  
% Sea: Season
% IndxSea{1} = [IndxMon{1};IndxMon{2};IndxMon{12}]; % Winter season (DJF)
  IndxSea{1} = [IndxMon{1};IndxMon{2};IndxMon{3}; ...
      IndxMon{10};IndxMon{11};IndxMon{12}]; % Winter season (Jing et al. 2019)
  
jultime_wt = jultime(IndxSea{1});
tau_llav_wt= tau_llav(IndxSea{1});
tau_wt     = tau_raw(:,:,IndxSea{1});
Ro_wt      = Ro_raw(:,:,IndxSea{1});
Ro_wt_std = nanstd(Ro_wt,0,3);
% Ro_wt_av   = nanmean(Ro_wt,3);
% kesai_wt   = kesai_raw(:,:,IndxSea{1});
% figure;hist(Ro_wt_av(indxRo))
% Ro_wt_av_Ro=Ro_wt_av(indxRo);
% sum(abs(Ro_wt_av_Ro)<0.15)./length(Ro_wt_av_Ro)
w_nl_wt = w_nl_raw(:,:,IndxSea{1});
w_li_wt = w_li_raw(:,:,IndxSea{1});
% w_df_wt = w_df_raw(:,:,IndxSea{1});
w_ce_wt = w_ce_raw(:,:,IndxSea{1});
t_ce_wt = t_ce_raw(:,:,IndxSea{1});  
  
tau_st_limt = mean(tau_llav_wt); 
tau_ns_limt = mean(tau_llav_wt);
% ---------------------------

% --- find winter storm days ---
indx_st = find(tau_llav_wt>=tau_st_limt); 
disp(['# of winter storm days: ',num2str(length(indx_st))])

tau_st     = tau_wt(:,:,indx_st);
tau_st_av  = mean(tau_st,3);
Ro_st      = Ro_wt(:,:,indx_st);
Ro_st_av   = nanmean(Ro_st,3);
w_li_st    = w_li_wt(:,:,indx_st);
w_li_st_av = nanmean(w_li_st,3);
w_nl_st    = w_nl_wt(:,:,indx_st);
w_nl_st_av = nanmean(w_nl_st,3);
% w_df_st = w_df_wt(:,:,indx_st);
% w_df_st_av = nanmean(w_df_st,3);
w_ce_st = w_ce_wt(:,:,indx_st);
w_ce_st_av = nanmean(w_ce_st,3);
% kesai_st = kesai_wt(:,:,indx_st);
% kesai_std_st = nanstd(kesai_st,0,3);
t_ce_st = t_ce_wt(:,:,indx_st);

[Q_li_st,~,~] = CalcVerticalEddyHeatFluxFunc(w_li_st,t_ce_st,rho_w);
[Q_nl_st,~,~] = CalcVerticalEddyHeatFluxFunc(w_nl_st,t_ce_st,rho_w);
[Q_ce_st,~,~] = CalcVerticalEddyHeatFluxFunc(w_ce_st,t_ce_st,rho_w);
% ----------------------

% --- find non winter storm days ---
indx_ns = find(tau_llav_wt<tau_ns_limt);
disp(['# of non winter storm days: ',num2str(length(indx_ns))])

tau_ns = tau_wt(:,:,indx_ns);
tau_ns_av = mean(tau_ns,3);
Ro_ns      = Ro_wt(:,:,indx_ns);
Ro_ns_av   = nanmean(Ro_ns,3);
w_li_ns = w_li_wt(:,:,indx_ns);
w_li_ns_av = nanmean(w_li_ns,3);
w_nl_ns = w_nl_wt(:,:,indx_ns);
w_nl_ns_av = nanmean(w_nl_ns,3);
% w_df_ns = w_df_wt(:,:,indx_ns);
% w_df_ns_av = nanmean(w_df_ns,3);
w_ce_ns = w_ce_wt(:,:,indx_ns);
w_ce_ns_av = nanmean(w_ce_ns,3);
t_ce_ns = t_ce_wt(:,:,indx_ns);

[Q_li_ns,~,~] = CalcVerticalEddyHeatFluxFunc(w_li_ns,t_ce_ns,rho_w);
[Q_nl_ns,~,~] = CalcVerticalEddyHeatFluxFunc(w_nl_ns,t_ce_ns,rho_w);
[Q_ce_ns,~,~] = CalcVerticalEddyHeatFluxFunc(w_ce_ns,t_ce_ns,rho_w);
% ---------------------------------
    
% --- compute the difference in W and Q from storm-free to storm days ---
w_ce_df = w_ce_st_av - w_ce_ns_av; % nanmean(Q_ce_df(:)) = 13.0692
w_nl_df = w_nl_st_av - w_nl_ns_av; % nanmean(Q_nl_df(:)) = 2.3855
w_li_df = w_li_st_av - w_li_ns_av; % nanmean(Q_li_df(:)) = 0.6762
w_ce_df_Ro = w_ce_df(Ro_wt_std > Ro_std_limit); %nanmean(w_ce_df_Ro)=46.1787 
w_nl_df_Ro = w_nl_df(Ro_wt_std > Ro_std_limit); %nanmean(w_nl_df_Ro)=20.0712
w_li_df_Ro = w_li_df(Ro_wt_std > Ro_std_limit); %nanmean(w_li_df_Ro)=8.2646


Q_ce_df = Q_ce_st - Q_ce_ns; % nanmean(Q_ce_df(:)) = 13.0692
Q_nl_df = Q_nl_st - Q_nl_ns; % nanmean(Q_nl_df(:)) = 2.3855
Q_li_df = Q_li_st - Q_li_ns; % nanmean(Q_li_df(:)) = 0.6762
Q_ce_df_Ro = Q_ce_df(Ro_wt_std > Ro_std_limit); %nanmean(Q_ce_df_Ro)=46.1787 
Q_nl_df_Ro = Q_nl_df(Ro_wt_std > Ro_std_limit); %nanmean(Q_nl_df_Ro)=20.0712
Q_li_df_Ro = Q_li_df(Ro_wt_std > Ro_std_limit); %nanmean(Q_li_df_Ro)=8.2646
  
% indxRo = find(Ro_wt_std>Ro_std_limit); % length(indxRo) = 8235
% w_li_Ro = w_li_st_av(indxRo);
% w_nl_Ro = w_nl_st_av(indxRo); 
% w_ce_Ro = w_ce_st_av(indxRo); 
% Q_li_Ro = Q_li_st_av(indxRo);
% Q_nl_Ro = Q_nl_st_av(indxRo); 
% Q_ce_Ro = Q_ce_st_av(indxRo);  
% --------------------------------------

% --- compute correlation ---
% despike W and Q
fillmethod = NaN;
findmethod = 'median';

[Q_li_df_Ro_Dspk,idx,~,~,~] = filloutliers(Q_li_df_Ro,fillmethod,findmethod);
[Q_nl_df_Ro_Dspk,idx,~,~,~] = filloutliers(Q_nl_df_Ro,fillmethod,findmethod);
[Q_ce_df_Ro_Dspk,idx,~,~,~] = filloutliers(Q_ce_df_Ro,fillmethod,findmethod);

nanstd(Q_ce_df_Ro_Dspk) % nanvar(Q_ce_df_Ro_Dspk) = 7.8061e+03
nanstd(Q_nl_df_Ro_Dspk) % nanvar(Q_nl_df_Ro_Dspk) = 1.2100e+03
nanvar(Q_li_df_Ro_Dspk) % nanvar(Q_li_df_Ro_Dspk) = 138.4845
nanvar((Q_nl_df_Ro_Dspk-nanmean(Q_ce_df_Ro_Dspk)))
nanvar((Q_nl_df_Ro_Dspk-Q_ce_df_Ro_Dspk))

figure;
plot(Q_nl_df_Ro_Dspk,'b');hold on;plot(Q_ce_df_Ro_Dspk,'r');grid on;
figure;
plot(Q_nl_df_Ro_Dspk,Q_ce_df_Ro_Dspk,'b*');hold on;plot(Q_ce_df_Ro_Dspk,'r');grid on;


% 
% sz = sine
% ind2sub();
% 
% figure
%  subplot(4,1,1);
%   plot(Q_nl_df_Ro,'b');hold on;plot(Q_nl_df_Ro_Dspk,'r');grid on;
%  subplot(4,1,2)
%   plot(Q_ce_df_Ro,'b');hold on;plot(Q_ce_df_Ro_Dspk,'r');grid on;
%  subplot(4,1,3)
%   plot(Q_nl_df_Ro_Dspk,'b');hold on;plot(Q_ce_df_Ro_Dspk,'r');grid on;
%  subplot(4,1,4)
%   plot(Q_nl_df_Ro_Dspk./Q_ce_df_Ro_Dspk,'b');hold on;plot(Q_ce_df_Ro_Dspk,'r');grid on;  
%   set(gca,'ylim',[-5 5],'ytick',[-5:1:5]);grid on;
%   
%   k = [500:700];
%   [rQ_nl2ce,P] = corrcoef(Q_nl_df_Ro_Dspk(k),Q_ce_df_Ro_Dspk(k));
%   
%   
% % Q_nl_df_Ro(abs(Q_nl_df_Ro-nanmean(Q_nl_df_Ro))>3.*nanstd(Q_nl_df_Ro))=nan;
% % Q_ce_df_Ro(abs(Q_ce_df_Ro-nanmean(Q_ce_df_Ro))>3.*nanstd(Q_ce_df_Ro))=nan;
% 
% % find No-NaN values
% indxNoNaN_Wnl = find(~isnan(w_nl_df_Ro));
% indxNoNaN_Wli = find(~isnan(w_li_df_Ro));
% indxNoNaN_Wce = find(~isnan(w_ce_df_Ro));
% % union([1 2 3], [2 3 4 5 6])
% indxNoNaN_W = intersect(intersect(indxNoNaN_Wnl,indxNoNaN_Wli),indxNoNaN_Wce);
% 
% Q_nl_df_Ro_Dspk(abs(Q_nl_df_Ro_Dspk)<10)=nan;
% Q_ce_df_Ro_Dspk(abs(Q_ce_df_Ro_Dspk)<10)=nan;
% 
indxNoNaN_Qnl = find(~isnan(Q_nl_df_Ro_Dspk));
% indxNoNaN_Qli = find(~isnan(Q_li_df_Ro_Dspk));
indxNoNaN_Qce = find(~isnan(Q_ce_df_Ro_Dspk));
% indxNoNaN_Q = intersect(intersect(indxNoNaN_Qnl,indxNoNaN_Qli),indxNoNaN_Qce);
indxNoNaN_Q = intersect(indxNoNaN_Qnl,indxNoNaN_Qce); length(indxNoNaN_Q)
% 
% % compute correlation and P value
% [rW_nl2li,P] = corrcoef(w_nl_df_Ro(indxNoNaN_W),w_li_df_Ro(indxNoNaN_W));
% [rW_nl2ce,P] = corrcoef(w_nl_df_Ro(indxNoNaN_W),w_ce_df_Ro(indxNoNaN_W));
% [rQ_nl2li,P] = corrcoef(Q_nl_df_Ro_Dspk(indxNoNaN_Q),Q_li_df_Ro_Dspk(indxNoNaN_Q));
% [rQ_nl2ce,P] = corrcoef(Q_nl_df_Ro_Dspk(indxNoNaN_Q),Q_ce_df_Ro_Dspk(indxNoNaN_Q));

% rW_nl2li = corr2(w_nl_df_Ro(indxNoNaN_W),w_li_df_Ro(indxNoNaN_W)); % rW_nl2li=0.5174
% rW_nl2ce = corr2(w_nl_df_Ro(indxNoNaN_W),w_ce_df_Ro(indxNoNaN_W)); % rW_nl2ce=0.1728
% rQ_nl2li = corr2(Q_nl_df_Ro(indxNoNaN_Q),Q_li_df_Ro(indxNoNaN_Q)); % rQ_nl2li=0.5223
% rQ_nl2ce = corr2(Q_nl_df_Ro(indxNoNaN_Q),Q_ce_df_Ro(indxNoNaN_Q)); % rQ_nl2ce=0.3341
% =====================


%% === make pics ===
% --- plot Wvel before, during and after storms ---
% plt_indx = 51;
% W_lim_nl = [-2 2]; W_nl_Ticks =[-2:1:2];
% 
% figure
%  subplot(4,3,1)
%   pcolor(Lon,Lat,w_nl_raw(:,:,plt_indx-1).*86400);shading interp;caxis(W_lim_nl);
%  subplot(4,3,2) 
%   pcolor(Lon,Lat,w_nl_raw(:,:,plt_indx).*86400);shading interp;caxis(W_lim_nl);
%  subplot(4,3,3)
%   pcolor(Lon,Lat,w_nl_raw(:,:,plt_indx+2).*86400);shading interp;caxis(W_lim_nl);
%  subplot(4,3,4)
%   pcolor(Lon,Lat,w_ce_raw(:,:,plt_indx-1).*86400);shading interp;caxis(W_lim_nl);
%  subplot(4,3,5) 
%   pcolor(Lon,Lat,w_ce_raw(:,:,plt_indx).*86400);shading interp;caxis(W_lim_nl);
%  subplot(4,3,6)
%   pcolor(Lon,Lat,w_ce_raw(:,:,plt_indx+2).*86400);shading interp;caxis(W_lim_nl);
%   
%   polarmap
% ------------------------

% --- compute ocean and atmosphere heat flux ---
% [Q_a_st,~,~] = CalcVerticalEddyHeatFluxFunc(w_nl_st,t_ce_st,rho_w);
% [Q_o_st,~,~] = CalcVerticalEddyHeatFluxFunc(w_ce_st-w_nl_st,t_ce_st,rho_w);
% [Q_t_st,~,~] = CalcVerticalEddyHeatFluxFunc(w_ce_st,t_ce_st,rho_w);
% 
% figure
% Q_lim_nl=[-200 200];
%  subplot(2,3,4)
%   pcolor(Lon,Lat,Q_a_st);shading interp;caxis(Q_lim_nl);colorbar;
%   title(['Qatm, av:',num2str(nanmean(Q_a_st(:))),' W/m^2']);ylabel('winter storms')
%  subplot(2,3,5) 
%   pcolor(Lon,Lat,Q_o_st);shading interp;caxis(Q_lim_nl);colorbar;
%   title(['Qocnv av:',num2str(nanmean(Q_o_st(:))),' W/m^2'])
%  subplot(2,3,6)
%   pcolor(Lon,Lat,Q_t_st);shading interp;caxis(Q_lim_nl);colorbar;
%   title(['Qocnv av:',num2str(nanmean(Q_t_st(:))),' W/m^2'])
%   polarmap
%   
%   nanmean(Q_t_st(:)), nanmean(Q_a_st(:)), nanmean(Q_o_st(:))
%   nanvar(Q_t_st(:)), nanvar(Q_a_st(:)), nanvar(Q_o_st(:)), nancov(Q_a_st(:),Q_o_st(:))
%   
% A = randn(10,1);
% B = randn(10,1);
% R = corrcoef(A+B,A+B)
% R = corrcoef(A,A+B)
% R = corrcoef(B,A+B)
% -------------------  

pltsaveTauWvelContoursSubScript % plot and save Tau, W contours

pltsaveQContoursSubScript % plot and save Q contours

% ======================


%% === output data ===
print(f2,'-dpng','-r500',pic2) 
print(f3,'-dpng','-r500',pic3) 
% ====================