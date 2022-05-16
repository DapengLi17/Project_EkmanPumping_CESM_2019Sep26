%% ===readme===

% descrip: matlab scripts plot near bottom water temperature,
% N2 (buoyancy frequency squared), S2 (vel shear squared), 
% Ri (gradient Richardson number), epsilon (dissipation rates) 
% Krho (turbulent diffusivity) and I (isotropy parameters)
% during 2017July KW.  

% update history:
% v1.0 DL 2019Feb14
% v1.1 DL 2019Mar25

% extra notes:
% near bottom(BBL) is ~ 1.5-6 MAB
% =============


%% ====set up environments====
clear all;close all;clc;

  date_str='2019Mar25';
  
inputfile1='../data_after_manipulation/RBR_Twater_15-28July2017KW_2018Sep07.mat';
inputfile2='../data_after_manipulation/Signature_burst-av_enu_vel_15-28July2017KW_2018Oct11.mat';  
inputfile3='../data_after_manipulation/Eps_Signature_HR_15-28July2017KW_2019Feb13.mat';

pic1=['../pics/TwN2S2RiEpsKrho_BBL_Contours_',date_str,'.png'];
pic2=['../pics/RiKrho_BBL_ScatterPlot_',date_str,'.png'];

% ~~~ study period ~~~  
% whole depl period, see log book    
   stime=datenum([2017 07 15 00 00 00]);
   ftime=datenum([2017 07 28 00 00 00]);

% picture period
   jultime_plt=[datenum([2017 07 15 00 00 00]) ...
       datenum([2017 07 28 00 00 00])];

% Spring/Neap tides time  
  Np_time=[datenum([2017 07 16 00 00 00]) datenum([2017 07 19 00 00 00])];
  Sp_time=[datenum([2017 07 23 00 00 00]) datenum([2017 07 26 00 00 00])];   
   
M_V=1.2e-6; % M_V=Molecular Viscosity of sea water  
lat  = 29.1717;  
dpth_av=22.9; % 22.9 m av wate depth 
av_method='bootstrap';
%=============================


%% ====load data==============
RBR=load(inputfile1);
Burst=load(inputfile2); % Signature Burst
HR=load(inputfile3);
%=============================


%% ===== data analyzation ===
% ============== 1 RBR Twater and N2 ================
 RBR_cell_indx=[6 8];

MAB_RBR=RBR.MAB_RBR(RBR_cell_indx);
Tw_raw=RBR.Tw_av(RBR_cell_indx,:);

disp('1/6 hr av')
[Tw_av,~,~,jultime_Tw_av]=multi_layers_time_av_btsp_mn(stime, ...
    ftime,Tw_raw,RBR.jultime_av,av_method);
 Tw_av_sort=sort(Tw_av,'descend');
 
Sal = 42.8.*ones(size(Tw_av_sort));
P=sw_pres(dpth_av-MAB_RBR,lat);
Press=repmat(P,1,size(Tw_av_sort,2));
[N2_raw,q,p_ave] = sw_bfrq(Sal,Tw_av_sort,Press,lat);

Tw_Diff=abs(diff(Tw_av_sort)); % hist(Tw_Diff)
% The RBRsolo³ T is calibrated to an accuracy of ±0.002ºC 
% https://rbr-global.com/products/compact-loggers/rbrsolo-t
% N2_raw(Tw_Diff<0.002)=nan;nanmin(N2_raw(:))=1.8989e-06;
% For Tw_Diff>=0.002, the min N2 resoloved is 1.8989e-06 s^(-2), 
% so I assign N2_raw(Tw_Diff<0.002) is 2e-6 s^(-2); 
N2_raw(Tw_Diff<0.002)=2e-6;

MAB_N2=dpth_av-sw_dpth(p_ave(:,1),lat);

disp('1 hr av')
[N2_av,~,~,jultime_N2_av]=multi_layers_time_av_btsp_mn(stime, ...
    ftime,N2_raw,jultime_Tw_av,av_method);
% =============================================

% ============== 2 compute S2 ==============
 bin_indx_Burst=3:12;
 MAB_S2=Burst.MAB_burst(bin_indx_Burst)'; % flipud([1;2;3;4;5])
 % MAB_S2([1 end])
 U1_brst_av=Burst.U1_brst_av(bin_indx_Burst,:);
 U2_brst_av=Burst.U2_brst_av(bin_indx_Burst,:);

% ~~~ use central difference method to compute horizontal vel shear ~~~
[dudz,dvdz] = calc_HorizVelShear(U1_brst_av,U2_brst_av,MAB_S2);
  S2_raw=dudz.^2+dvdz.^2;
% min S2 resolved: 1cm/s(0.01m/s) vel spd diff over 1 cell(0.5m)
% (0.01/0.5).^2 = 4e-4
% S2_raw(S2_raw<1e-4)=nan;
S2_raw(S2_raw<1e-4)=1e-4;

disp('1 hr av')
[S2_av,~,~,jultime_S2_av]=multi_layers_time_av_btsp_mn(stime, ...
    ftime,S2_raw,Burst.jultime_brst_av,av_method);

disp('1 hr av')
[S2_dpth_time_av,~,~,~]=multi_layers_time_av_btsp_mn(stime, ...
    ftime,nanmean(S2_raw),Burst.jultime_brst_av,av_method);
% ====================================

% =========== 3 compute Ri ===========
% MAB_S2
disp('1/6 hr av')
[S24Ri_raw,~,~,~]=multi_layers_time_av_btsp_mn(stime, ...
   ftime,nanmean(S2_raw),Burst.jultime_brst_av,av_method);

Ri_raw=N2_raw./S24Ri_raw;

disp('1hr av')
[Ri_av,~,~,jultime_Ri_av]=multi_layers_time_av_btsp_mn(stime, ...
    ftime,Ri_raw,jultime_Tw_av,av_method);
% ===================================

% =========== 4 Signature HR Epsilon =============
% bin_indx_HR=1:size(HR.eps_raw,1);
bin_indx_HR=1:100;
eps_raw=HR.eps_raw(bin_indx_HR,:);
MAB_HR=HR.MAB_HR(bin_indx_HR);

% vertical bin av, 4bin = 4 * 0.05 = 0.2m
for i = 1 : size(eps_raw,2); % loop over time
   eps_dpth_av(:,i)=binave_aa(eps_raw(:,i),4); 
end
% binave_aa([1 2 nan 4], 2)
MAB_eps=binave_aa(MAB_HR,4);

% figure;
%  [c,h]=contourf(HR.jultime_HR,MAB_eps,log10(eps_dpth_av));
%   set(h,'linestyle','none');caxis([-9 -5]);
%   datetick('x','dd/mmm')

disp('1 hr av')
 [eps_time_dpth_av,~,~,jultime_eps_av]= ...
     multi_layers_time_av_btsp_mn(stime, ...
    ftime,eps_dpth_av,HR.jultime_HR,av_method);

% I found and removed a spike of epsilon 
% between 2017 July 15 @ 22:00 to 16 @ 02:00
IndxClean=find(jultime_eps_av>=datenum([2017 07 15 22 00 00]) & ...
    jultime_eps_av<=datenum([2017 07 16 02 00 00]));
[m,n]=find(eps_time_dpth_av(:,IndxClean)>1e-5);
eps_time_dpth_av(m,IndxClean(n))=nan;

Eps_DpthAv=nan(1,size(eps_time_dpth_av,2));
for i = 1 : size(eps_time_dpth_av,2) % loop over time
 if sum(~isnan(eps_time_dpth_av(:,i)))>9; % eps_time_dpth_av    25x312    
     % at least average over 10 points
     [~,Eps_DpthAv(i),~] = bootstrap5(eps_time_dpth_av(:,i));
 end
end
   
disp('24 hr av')
 [Eps_DpthTimeAv,Eps_DpthTimeMin,Eps_DpthTimeMax,jultime_Eps_DpthTimeAv] = ...
     multi_layers_time_av_btsp_mn(stime,ftime,Eps_DpthAv, ...
     jultime_eps_av,av_method);
% ================================================

% ===== 5 compute I and Krho =====
 I=Eps_DpthAv./M_V./N2_av; % Isotropy parameters, see Thorpe 2007
 % hist(log10(I))
 
 % sum(I>200)./length(I) = 0.9391
 % sum(I>100)./length(I) = 0.9583
 
 Krho=2.*M_V.*I.^(0.5); % diapycnal diffusivity in energetic regime, see
% Shil et al. 2005 Parametrization of turbulent fluxes and scales using
% homogeneous sheared stably stratified turbulence simulations, Table 3.  
 
% --- error analysis on Osborn 1980 model --- 
  Krho_OS80=0.2.*Eps_DpthAv./N2_av; % Osborn 1980 model
% plot(jultime_N2_av,Krho,'b');hold on;grid on;
% plot(jultime_N2_av,Krho_OS80,'r');set(gca,'yscale','log');tlabel;
[ratio_min,ratio_mn,ratio_max]=bootstrap5(Krho_OS80./Krho); 
% ratio_min=19.1196; ratio_mn=20.9968; ratio_max = 22.9396;0.
 
% --- correlation coefficient between Ri and Krho ---
% IndxNoNaN=find(~isnan(Krho));
% corrcoef(log10(Ri_av(IndxNoNaN)),log10(Krho(IndxNoNaN)))
%     1.0000   -0.8705
%    -0.8705    1.0000
% ---------------------------------------------------

% --- fit Krho with Ri ---
% Krho=a.*Ri.^(b), format follow Peters et al. 1988 On the parametrization 
% of equatoaial turbulence, see Eq (10)-(11)
  [B,stats] =robustfit(log10(Ri_av(:)),log10(Krho(:)));
    % stats.coeffcorr = 0.9661, correlation of fit is 97%
    a=10.^(B(1)); % a=9.4542e-06; 
    b=B(2); % b= -0.9769
  Ri_fit=[4e-3 1e-2 1e-2 4e-1];
  Krho_fit=a.*Ri_fit.^b; 
  Krho_fit_min=a.*Ri_fit.^(b-stats.se(2)); 
  Krho_fit_max=a.*Ri_fit.^(b+stats.se(2)); 
 
% Peters et al. 1988 parametrizations  
Ri_Peters=[0.25 0.4];
%   Kh_Peters_av=3.0e-9.*Ri_Peters.^(-9.6);
%   Kh_Peters_min=3.0e-9.*Ri_Peters.^(-9.6-3.7);
%   Kh_Peters_max=3.0e-9.*Ri_Peters.^(-9.6+3.7);
Krho_Peters=1.1e-8.*Ri_Peters.^(-9.2);
% ==============================================


% ===== 5 find the strong daily mixing between 22–26July2017 =====
Indx_Ebbs=find(jultime_eps_av>=datenum([2017 07 22 00 00 00]) & ...
    jultime_eps_av<=datenum([2017 07 26 00 00 00]));
jultime_Ebbs=jultime_eps_av(Indx_Ebbs); % datestr(jultime_Ebbs)
Eps_Ebbs=Eps_DpthAv(Indx_Ebbs);

Eps_reshp=reshape(Eps_Ebbs,24,numel(Eps_Ebbs)./24);
jultime_Ebbs_reshp=reshape(jultime_Ebbs,24,numel(Eps_Ebbs)./24);
[~,IndxEpsMax]=max(Eps_reshp);

for i = 1 : length(IndxEpsMax)
  jultime_EpsMax(i)=jultime_Ebbs_reshp(IndxEpsMax(i),i);
end
% % ===========================================


%% === make pics ===
% contour plot
f1=figure('Renderer','painters');

  % set pic size
  pic1_size=[6.62 7]; % pic size: unit [inch] 3.155 inch is the width(length in x direction) of JPO 1 column pic, height varies according to your needs
  set(f1,'Units','inches','Position',[5,5,pic1_size]);
  font_size=9;
  ylabel_position=[-0.07, 0.5, 0];
  Xlim=jultime_plt;% datevec(Xlim)
  xticks=linspace(Xlim(1),Xlim(2),14); % datevec(xticks)  time interval 
 
% ~~~ generate subplot position ~~~  
  row_num=6;col_num=1;margin_left=0.1;
  margin_right=0.1;margin_top=0.015;margin_botm=0.07;
  pics_dist_x=0;pics_dist_y=0.03;
 
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

subplot('position',sbplt_posit(1,:)); % Tw
  h1=plot(jultime_Tw_av,Tw_av_sort(1,:),'r','linewidth',1.5);hold on;
  h2=plot(jultime_Tw_av,Tw_av_sort(2,:),'b','linewidth',1.5);grid on;
    set(gca,'xlim',Xlim,'XTick',xticks,'XTickLabel',[],'ylim',[32 34.6], ...
      'ytick',[32:1:34],'fontsize',font_size);
    ylabel('T [\circC]','Units','Normalized', ...
        'Position',ylabel_position,'fontsize',font_size);
    hold on; % plot spring/neap tide time
  plot(Sp_time,[34.4 34.4],'r','linewidth',2);
  plot(Np_time,[34.4 34.4],'b','linewidth',2);
    legend([h1,h2],{'8 MAB','0.72 MAB'},'orientation','horizontal', ...
    'Position',[0.61 0.87 0.27 0.03],'fontsize',font_size-1);legend boxoff;
    text(jultime_N2_av(5),34.3,'a','Fontsize',font_size);
         
subplot('position',sbplt_posit(2,:)); % depth-av N2, S2     
  [Ax,h1,h2]=plotyy(jultime_N2_av,N2_av,jultime_S2_av,S2_dpth_time_av);
     set(h1,'linewidth',1.2,'color','b');
     set(h2,'linewidth',1.2,'color','r');
     set(get(Ax(1),'Ylabel'),'String','N^2 [s^{-2}]','FontSize',font_size);
     set(get(Ax(1),'Ylabel'),'Units', 'Normalized', 'Position', ylabel_position);
     set(get(Ax(2),'Ylabel'),'String','S^2 [s^{-2}]','FontSize',font_size);    
     set(Ax(1),'xlim',Xlim,'XTick',xticks,'XTickLabel',[],'yscale','log', ...
       'ylim',[1e-6 1.2e-3],'ytick',[1e-6 1e-5 1e-4 1e-3], ...
       'ycolor','b','fontsize',font_size-1);
     set(Ax(2),'xlim',Xlim,'XTick',xticks,'XTickLabel',[],'yscale','log', ...
       'ylim',[1e-6 1.2e-3],'ytick',[1e-6 1e-5 1e-4 1e-3], ...
       'ycolor','r','fontsize',font_size-1);grid on;  
     text(jultime_N2_av(5),6e-4,'b','Fontsize',font_size);

subplot('position',sbplt_posit(3,:)); % Ri and Krho
 [Ax,h1,h2]=plotyy(jultime_Ri_av,Ri_av,jultime_N2_av,Krho);grid on;
   set(h1,'linewidth',1.2,'color','b');
   set(h2,'linewidth',1.2,'color','r');
%    set(Ax(1),'yscale','log','Ylim',[1e-3 1e0], ...
%        'ytick',[1e-3 1e-2 1e-1 1e0],'Fontsize',font_size-1);
%    set(Ax(2),'yscale','log','Ylim',[1e-5 1e-2], ...
%        'ytick',[1e-5 1e-4 1e-3 1e-2],'Fontsize',font_size-1);
   set(Ax(1),'yscale','log','Ylim',[1e-5 1e0], ...
       'ytick',[1e-4 1e-2 1e0],'Fontsize',font_size-1);
   set(Ax(2),'yscale','log','Ylim',[1e-5 1e0], ...
       'ytick',[1e-4 1e-2 1e0],'Fontsize',font_size-1);
   set(get(Ax(1),'Ylabel'),'String','Ri','FontSize',font_size);
   set(get(Ax(1),'Ylabel'),'Units','Normalized', ...
       'Position',ylabel_position);
   set(get(Ax(2),'Ylabel'),'String','K_{\rho} [m^2/s]', ...
       'FontSize',font_size);
   set(Ax(1),'Xlim',Xlim,'xtick',xticks,'xticklabel',[],'ycolor','b');
   set(Ax(2),'Xlim',Xlim,'xtick',xticks,'xticklabel',[],'ycolor','r');
%   plot(jultime_plt,[1 1],'r','linewidth',1.5);
   text(jultime_N2_av(5),5.e-1,'c','Fontsize',font_size); 
    
subplot('position',sbplt_posit(4,:)); % S2 contours
% min(log10(S2_av(:))),max(log10(S2_av(:)))
 [c,h]=contourf(jultime_S2_av,MAB_S2,log10(S2_av));
  set(h,'linestyle','none');caxis([-4 -3]);
%   hold on;
%  [c,h]=contour(jultime_S2_av,MAB_S2,log10(S2_av), ...
%      [-3.2,-3.2],'w');  
   set(gca,'xlim',Xlim,'XTick',xticks,'XTickLabel',[], ...
       'ylim',[1.5 7],'ytick',[2:2:6],'fontsize',font_size);
   ylabel(['log_{10}(S^2)',sprintf('\n'),'MAB'],'Units','Normalized', ...
        'Position',[-0.04 0.5 0],'fontsize',font_size);    
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.94 sbplt_posit(4,2) .01 sbplt_posit(4,4)-0.04], ...
       'ylim',[-4 -3],'YTick',[-4:.5:-3]);
   title(h1,'[s^{-2}]','FontSize',font_size-2);%colormap(sbplt3,jet)
   text(jultime_N2_av(5),6.4,'d','fontsize',font_size); 
   
subplot('position',sbplt_posit(5,:)); % Epsilon contours
% min(log10(eps_time_dpth_av(:))),max(log10(eps_time_dpth_av(:))),
 [c,h]=contourf(jultime_eps_av,MAB_eps,log10(eps_time_dpth_av));
   set(h,'linestyle','none');caxis([-9 -5]);hold on;
 [c,h]=contour(jultime_eps_av,MAB_eps,log10(eps_time_dpth_av), ...
     [-6,-6],'r');  
 plot(jultime_EpsMax,7.*ones(size(jultime_EpsMax)),'kv','MarkerFaceColor','k'); 
   set(gca,'xlim',Xlim,'XTick',xticks,'XTickLabel',[], ...
     'ylim',[1.5 7],'ytick',[2:2:6],'FontSize',font_size); 
   ylabel(['log_{10}(\epsilon)',sprintf('\n'),'MAB'],'Units','Normalized', ...
     'Position',[-0.04 0.5 0],'fontsize',font_size);
   
% add colorbar 
   h1=colorbar('FontSize',font_size-2);
   set(h1,'position',[0.94 sbplt_posit(5,2) .01 sbplt_posit(5,4)-0.04], ...
       'ylim',[-9 -5],'YTick',[-9:2:-5]);
   title(h1,'[W/kg]','FontSize',font_size-2);  
   text(jultime_N2_av(5),6.6,'e','fontsize',font_size); 

subplot('position',sbplt_posit(6,:)); % Depth-av Epsilon   
 h1=plot(jultime_eps_av,Eps_DpthAv,'b-','linewidth',1.2);hold on;
 h2=plot(jultime_Eps_DpthTimeAv,Eps_DpthTimeAv,'ro-','linewidth',1.2, ...
     'MarkerSize',4);grid on; 
 h_fil=fill([jultime_Eps_DpthTimeAv(1) jultime_Eps_DpthTimeAv  ...
     fliplr(jultime_Eps_DpthTimeAv)],[Eps_DpthTimeMax(1)  ...
     Eps_DpthTimeMin fliplr(Eps_DpthTimeMax)],'r');
%  plot(jultime_EpsMax,3e-6.*ones(size(jultime_EpsMax)),'kv','MarkerFaceColor','k');
   set(h_fil,'FaceColor','r','EdgeAlpha',0.0,'FaceAlpha',0.2);
    legend([h1,h2],{'hourly av','daily av'},'FontSize',font_size-1, ...
     'orientation','horizontal','Position',[0.61 0.075 0.25 0.02]);legend boxoff;
   set(gca,'xlim',Xlim,'XTick',xticks,'yscale','log', ...
     'ylim',[1.5e-9 3e-6],'ytick',[1e-8 1e-7 1e-6], ...
     'FontSize',font_size-1); 
   ylabel(['\epsilon [W/kg]'],'Units','Normalized', ...
     'Position',ylabel_position,'fontsize',font_size);
   xlabel('Day/July/2017','FontSize',font_size); 
   text(jultime_N2_av(5),1.5e-6,'f','fontsize',font_size); 
   datetick('x','DD','keepticks','keeplimits'); 
% ======================


%% =============== scatter plot =================
f2=figure('Renderer','painters'); 

  % set pic size
%   pic2_size=[6.62 3];
  pic2_size=[3.155 3.155]; % pic size: unit [inch] 3.155 inch is the width(length in x direction) of JPO 1 column pic, height varies according to your needs
  set(f2,'Units','inches','Position',[5,5,pic2_size]);
  font_size=10;
  ylabel_position=[-0.15, 0.5, 0];
  
% ~~~ generate subplot position ~~~  
  row_num=1;col_num=1;margin_left=0.19;
  margin_right=0.05;margin_top=0.05;margin_botm=0.19;
  pics_dist_x=0.0;pics_dist_y=0.;
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
 subplot('position',sbplt_posit(1,:)); % Krho vs Ri
  h1=plot(Ri_av,Krho,'b*','markersize',3);hold on;
  h2=plot(Ri_fit,Krho_fit,'r-','linewidth',1.5);grid on;
  h3=plot(Ri_Peters,Krho_Peters,'k-','linewidth',1.5);% Peters et al. 1988 parametrization
  h_fil=fill([Ri_fit(1) Ri_fit fliplr(Ri_fit)], ...
    [Krho_fit_min(1) Krho_fit_min fliplr(Krho_fit_max)],'r');
  legend([h1,h2,h3],{'Obs.','Fit','Peters et al. 1988'}, ...
      'location','southwest','FontSize',font_size-2);legend boxoff;
    set(h_fil,'FaceColor','r','EdgeAlpha',0.0,'FaceAlpha',0.2);
    set(gca,'xscale','log','xlim',[1e-3 1e0],'xtick',[1e-3 1e-2 1e-1 1e0], ...
      'yscale','log','ylim',[1e-5 1e-2],'ytick',[1e-5 1e-4 1e-3 1e-2], ...
      'Fontsize',font_size);
%     text(1.3e-3,1.4e-5,'a','fontsize',font_size);
    xlabel('Ri','FontSize',font_size);         
    ylabel('K_{\rho} [m^2/s]','Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);

% % Peters et al. 1988 parametrization vs ours
% subplot('position',sbplt_posit(2,:)); hold on; box on;
% %  h1=plot(Ri_Peters,Kh_Peters_av,'b-','linewidth',1.5);
% %   h_fil=fill([Ri_Peters(1) Ri_Peters fliplr(Ri_Peters)], ...
% %       [Kh_Peters_min(1) Kh_Peters_min fliplr(Kh_Peters_max)],'b');
% %   set(h_fil,'FaceColor','b','EdgeAlpha',0.0,'FaceAlpha',0.2);
%  h2=plot(Ri_fit,Krho_fit,'r-','linewidth',1.5);
%   h_fil=fill([Ri_fit(1) Ri_fit fliplr(Ri_fit)], ...
%       [Krho_fit_min(1) Krho_fit_min fliplr(Krho_fit_max)],'r');
%   set(h_fil,'FaceColor','r','EdgeAlpha',0.0,'FaceAlpha',0.2);
%   legend([h2,h1],{'NAG 2017','Peters et al. 1988'}, ...
%       'location','northwest','FontSize',font_size);legend boxoff;
%   set(gca,'xscale','log','xlim',[1e-3 1e0],'xtick',[1e-3 1e-2 1e-1 1e0], ...
%       'yscale','log','ylim',[5e-7 5e-1],'ytick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1], ...
%       'Fontsize',font_size);grid on;
%   text(1.3e-3,1.e-6,'b','fontsize',font_size);
%   xlabel('Ri','FontSize',font_size);       

%  subplot('position',sbplt_posit(1,:)); % N2 vs Eps
%   plot(N2_av,Eps_DpthTimeAv,'bo','markersize',3);grid on;
%     set(gca,'xscale','log','yscale','log','xlim',[1e-6 5e-4], ...
%       'xtick',[1e-6 1e-5 1e-4],'ylim',[1e-9 1e-5], ...
%       'ytick',[1e-9 1e-8 1e-7 1e-6 1e-5],'FontSize',font_size);
%     ylabel('\epsilon [W/kg]','Units','Normalized', ...
%         'Position',ylabel_position,'fontsize',font_size); 
%     xlabel('N^2 [s^{-2}]','fontsize',font_size);  
%     text(1.2e-6,2e-9,'a','fontsize',font_size); 
%      
%  subplot('position',sbplt_posit(2,:)); % S2 vs Eps
%   plot(S2_dpth_time_av,Eps_DpthTimeAv,'bo','markersize',3);grid on;
%    set(gca,'xscale','log','yscale','log','xlim',[7e-5 1.5e-3], ...
%        'xtick',[1e-4 1e-3],'ylim',[1e-9 1e-5], ...
%        'ytick',[1e-9 1e-8 1e-7 1e-6 1e-5],'yticklabel',[], ...
%        'FontSize',font_size);
%    xlabel('S^2 [s^{-2}]','fontsize',font_size);  
%    text(0.8e-4,2e-9,'b','fontsize',font_size); 
%     
%  subplot('position',sbplt_posit(3,:)); % Ri vs Eps
%   plot(Ri_av,Eps_DpthTimeAv,'bo','markersize',3);grid on;
%    set(gca,'xscale','log','yscale','log','xlim',[1e-3 1e0], ...
%        'xtick',[1e-3 1e-2 1e-1 1e0],'ylim',[1e-9 1e-5], ...
%        'ytick',[1e-9 1e-8 1e-7 1e-6 1e-5],'yticklabel',[], ...
%        'FontSize',font_size);
%    xlabel('Ri','fontsize',font_size);  
%    text(1.3e-3,2e-9,'c','fontsize',font_size);
% ==================================


%% === output files ===
Image_out=fillPage(f1,'margins',[0 0 0 0],'papersize',pic1_size,'inches');
% no margin
 print(f1,'-dpng','-r500',pic1)
 set(f1,Image_out); % do not edit the picture while matlab saves it
 
Image_out=fillPage(f2,'margins',[0 0 0 0],'papersize',pic2_size,'inches');
% no margin
 print(f2,'-dpng','-r500',pic2)
 set(f2,Image_out); % do not edit the picture while matlab saves it 
% =====================      