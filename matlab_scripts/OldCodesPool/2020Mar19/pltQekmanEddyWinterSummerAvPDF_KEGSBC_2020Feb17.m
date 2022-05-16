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

   date_str='2020Feb17';
   
infile1 = '../data_after_manipulation/QekmanEddyGlobalBorealSummerAv_2020Feb11.mat';
infile2 = '../data_after_manipulation/QekmanEddyGlobalBorealWinterAv_2020Feb11.mat';

pic1 = ['../pics/QekmanEddyWinterSummerAvPDF_KEGSBC_',date_str,'.png'];

% Constants4CESM_Global
addpath(genpath('Func4EkmanProject/'))

% Kuroshio Extension (KE)
   latKE_limits = [25 45]; 
   lonKE_limits = [130 170];
% Gulf Stream (GS)   
   latGS_limits = [25 50];  
   lonGS_limits = [275 325];
% Brazil Current (BC)   
   latBC_limits = [-55 -30];  
   lonBC_limits = [290 320];
   
QbinsCenter = [-10:10:40]; % for PDF
%=============================


%% === load data ===
Su = load(infile1);
Wt = load(infile2);
% ==================


%% === data analysis ===
lat_1d_raw = Su.lat_1d;
lon_1d_raw = Su.lon_1d;

% --- Kuroshio Extension (KE) ---
indxLatKE = find(lat_1d_raw >= latKE_limits(1) & lat_1d_raw <= latKE_limits(2));
lat_1dKE = lat_1d_raw(indxLatKE);
indxLonKE = find(lon_1d_raw >= lonKE_limits(1) & lon_1d_raw <= lonKE_limits(2));
lon_1dKE = lon_1d_raw(indxLonKE);

Q_li_wtKE = Wt.Q_li(indxLatKE,indxLonKE); 
Q_li_suKE = Su.Q_li(indxLatKE,indxLonKE); 
Q_nl_wtKE = Wt.Q_nl(indxLatKE,indxLonKE);
Q_nl_suKE = Su.Q_nl(indxLatKE,indxLonKE);
Q_ce_wtKE = Wt.Q_ce(indxLatKE,indxLonKE);
Q_ce_suKE = Su.Q_ce(indxLatKE,indxLonKE);

[N_Q_li_suKE,C_Q_li_suKE] = hist(Q_li_suKE(:),QbinsCenter); 
P_Q_li_suKE = N_Q_li_suKE./sum(N_Q_li_suKE);

[N_Q_li_wtKE,C_Q_li_wtKE] = hist(Q_li_wtKE(:),QbinsCenter); 
P_Q_li_wtKE = N_Q_li_wtKE./sum(N_Q_li_wtKE);

[N_Q_nl_suKE,C_Q_nl_suKE] = hist(Q_nl_suKE(:),QbinsCenter); 
P_Q_nl_suKE = N_Q_nl_suKE./sum(N_Q_nl_suKE);

[N_Q_nl_wtKE,C_Q_nl_wtKE] = hist(Q_nl_wtKE(:),QbinsCenter); 
P_Q_nl_wtKE = N_Q_nl_wtKE./sum(N_Q_nl_wtKE);

[N_Q_ce_suKE,C_Q_ce_suKE] = hist(Q_ce_suKE(:),QbinsCenter); 
P_Q_ce_suKE = N_Q_ce_suKE./sum(N_Q_ce_suKE);

[N_Q_ce_wtKE,C_Q_ce_wtKE] = hist(Q_ce_wtKE(:),QbinsCenter); 
P_Q_ce_wtKE = N_Q_ce_wtKE./sum(N_Q_ce_wtKE);
% ----------------------------------

% --- Gulf Stream (GS) ---
indxLatGS = find(lat_1d_raw >= latGS_limits(1) & lat_1d_raw <= latGS_limits(2));
lat_1dGS = lat_1d_raw(indxLatGS);
indxLonGS = find(lon_1d_raw >= lonGS_limits(1) & lon_1d_raw <= lonGS_limits(2));
lon_1dGS = lon_1d_raw(indxLonGS);

Q_li_wtGS = Wt.Q_li(indxLatGS,indxLonGS); 
Q_li_suGS = Su.Q_li(indxLatGS,indxLonGS); 
Q_nl_wtGS = Wt.Q_nl(indxLatGS,indxLonGS);
Q_nl_suGS = Su.Q_nl(indxLatGS,indxLonGS);
Q_ce_wtGS = Wt.Q_ce(indxLatGS,indxLonGS);
Q_ce_suGS = Su.Q_ce(indxLatGS,indxLonGS);

[N_Q_li_suGS,C_Q_li_suGS] = hist(Q_li_suGS(:),QbinsCenter); 
P_Q_li_suGS = N_Q_li_suGS./sum(N_Q_li_suGS);

[N_Q_li_wtGS,C_Q_li_wtGS] = hist(Q_li_wtGS(:),QbinsCenter); 
P_Q_li_wtGS = N_Q_li_wtGS./sum(N_Q_li_wtGS);

[N_Q_nl_suGS,C_Q_nl_suGS] = hist(Q_nl_suGS(:),QbinsCenter); 
P_Q_nl_suGS = N_Q_nl_suGS./sum(N_Q_nl_suGS);

[N_Q_nl_wtGS,C_Q_nl_wtGS] = hist(Q_nl_wtGS(:),QbinsCenter); 
P_Q_nl_wtGS = N_Q_nl_wtGS./sum(N_Q_nl_wtGS);

[N_Q_ce_suGS,C_Q_ce_suGS] = hist(Q_ce_suGS(:),QbinsCenter); 
P_Q_ce_suGS = N_Q_ce_suGS./sum(N_Q_ce_suGS);

[N_Q_ce_wtGS,C_Q_ce_wtGS] = hist(Q_ce_wtGS(:),QbinsCenter); 
P_Q_ce_wtGS = N_Q_ce_wtGS./sum(N_Q_ce_wtGS);
% -------------------------

% --- Brazil Current (BC) ---
indxLatBC = find(lat_1d_raw >= latBC_limits(1) & lat_1d_raw <= latBC_limits(2));
lat_1dBC = lat_1d_raw(indxLatBC);
indxLonBC = find(lon_1d_raw >= lonBC_limits(1) & lon_1d_raw <= lonBC_limits(2));
lon_1dBC = lon_1d_raw(indxLonBC);

Q_li_wtBC = Su.Q_li(indxLatBC,indxLonBC); 
Q_li_suBC = Wt.Q_li(indxLatBC,indxLonBC); 
Q_nl_wtBC = Su.Q_nl(indxLatBC,indxLonBC);
Q_nl_suBC = Wt.Q_nl(indxLatBC,indxLonBC);
Q_ce_wtBC = Su.Q_ce(indxLatBC,indxLonBC);
Q_ce_suBC = Wt.Q_ce(indxLatBC,indxLonBC);

[N_Q_li_suBC,C_Q_li_suBC] = hist(Q_li_suBC(:),QbinsCenter); 
P_Q_li_suBC = N_Q_li_suBC./sum(N_Q_li_suBC);

[N_Q_li_wtBC,C_Q_li_wtBC] = hist(Q_li_wtBC(:),QbinsCenter); 
P_Q_li_wtBC = N_Q_li_wtBC./sum(N_Q_li_wtBC);

[N_Q_nl_suBC,C_Q_nl_suBC] = hist(Q_nl_suBC(:),QbinsCenter); 
P_Q_nl_suBC = N_Q_nl_suBC./sum(N_Q_nl_suBC);

[N_Q_nl_wtBC,C_Q_nl_wtBC] = hist(Q_nl_wtBC(:),QbinsCenter); 
P_Q_nl_wtBC = N_Q_nl_wtBC./sum(N_Q_nl_wtBC);

[N_Q_ce_suBC,C_Q_ce_suBC] = hist(Q_ce_suBC(:),QbinsCenter); 
P_Q_ce_suBC = N_Q_ce_suBC./sum(N_Q_ce_suBC);

[N_Q_ce_wtBC,C_Q_ce_wtBC] = hist(Q_ce_wtBC(:),QbinsCenter); 
P_Q_ce_wtBC = N_Q_ce_wtBC./sum(N_Q_ce_wtBC);
% ----------------------------

clear Wt Su
% ======================


%% === make pics ===
f1=figure('Renderer','painters'); % Q

  pic1_size = [8.5 6]; 
  set(f1,'units','inches','position',[5,5,pic1_size])
%   set(f1,'units','normalized','position',[0,0,1,1])
  font_size = 10;
  y_limits = [0 0.8];
  y_ticks =  [0:0.2:0.8];
  x4text = -12;
  y4text = 0.75;
  ylabel_position=[-0.16, 0.5, 0];
  centerQ = C_Q_ce_suBC;
  
 % ~~~ generate subplot position ~~~  
  row_num=2;        col_num=3;
  margin_left=0.1;  margin_right=0.03;
  margin_top=0.05;  margin_botm=0.1;
  pics_dist_x=0.05; pics_dist_y=0.08;
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

   
  subplot('Position',sbplt_posit(1,:))
   hb = bar(centerQ,[P_Q_li_suKE;P_Q_nl_suKE;P_Q_ce_suKE]');grid on;
    set(hb(1),'FaceColor','b');set(hb(2),'FaceColor','r');
    set(hb(3),'FaceColor','k');
    legend(hb,{'Q_f','Q_{\xi}','Q_{eddy}^{50m}'},'location','northeast' ...
        );legend boxoff;
    set(gca,'fontsize',font_size,'ylim',y_limits,'ytick',y_ticks);
    title(['Kuroshio Extension'],'fontsize',font_size);
    ylabel(['Summer',sprintf('\n'),'[%]'],'Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    text(x4text,y4text,'a','fontsize',font_size,'FontWeight','bold');
   
  subplot('Position',sbplt_posit(2,:))
   hb = bar(centerQ,[P_Q_li_suGS;P_Q_nl_suGS;P_Q_ce_suGS]');grid on;
    set(hb(1),'FaceColor','b');set(hb(2),'FaceColor','r');
    set(hb(3),'FaceColor','k');
    set(gca,'fontsize',font_size,'ylim',y_limits,'ytick',y_ticks,'yticklabel',[]);
    title(['Gulf Stream'],'fontsize',font_size);
    text(x4text,y4text,'b','fontsize',font_size,'FontWeight','bold'); 

      
  subplot('Position',sbplt_posit(3,:))
   hb = bar(centerQ,[P_Q_li_suBC;P_Q_nl_suBC;P_Q_ce_suBC]');grid on;
    set(hb(1),'FaceColor','b');set(hb(2),'FaceColor','r');
    set(hb(3),'FaceColor','k');
    set(gca,'fontsize',font_size,'ylim',y_limits,'ytick',y_ticks,'yticklabel',[]);
    title(['Brazil Currents'],'fontsize',font_size);
    text(x4text,y4text,'c','fontsize',font_size,'FontWeight','bold'); 
    
  subplot('Position',sbplt_posit(4,:))
   hb = bar(centerQ,[P_Q_li_wtKE;P_Q_nl_wtKE;P_Q_ce_wtKE]');grid on;
    set(hb(1),'FaceColor','b');set(hb(2),'FaceColor','r');
    set(hb(3),'FaceColor','k');
    set(gca,'fontsize',font_size,'ylim',y_limits,'ytick',y_ticks);
    ylabel(['Winter',sprintf('\n'),'[%]'],'Units','Normalized', ...
      'Position',ylabel_position,'fontsize',font_size);
    xlabel(['[W/m^2]'],'fontsize',font_size);
    text(x4text,y4text,'d','fontsize',font_size,'FontWeight','bold');


  subplot('Position',sbplt_posit(5,:))
   hb = bar(centerQ,[P_Q_li_wtGS;P_Q_nl_wtGS;P_Q_ce_wtGS]');grid on;
    set(hb(1),'FaceColor','b');set(hb(2),'FaceColor','r');
    set(hb(3),'FaceColor','k');
    set(gca,'fontsize',font_size,'ylim',y_limits,'ytick',y_ticks,'yticklabel',[]);
    text(x4text,y4text,'e','fontsize',font_size,'FontWeight','bold');
    xlabel(['[W/m^2]'],'fontsize',font_size);
    
  subplot('Position',sbplt_posit(6,:))
   hb = bar(centerQ,[P_Q_li_wtBC;P_Q_nl_wtBC;P_Q_ce_wtBC]');grid on;
    set(hb(1),'FaceColor','b');set(hb(2),'FaceColor','r');
    set(hb(3),'FaceColor','k');
    set(gca,'fontsize',font_size,'ylim',y_limits,'ytick',y_ticks,'yticklabel',[]);
    text(x4text,y4text,'f','fontsize',font_size,'FontWeight','bold'); 
    xlabel(['[W/m^2]'],'fontsize',font_size);
% ====================

  
%% === output data ===  
export_fig(f1,pic1,'-r200','-nocrop')
% ====================
