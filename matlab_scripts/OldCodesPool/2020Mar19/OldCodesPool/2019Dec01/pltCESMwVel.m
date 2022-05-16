%% ===readme===

% descrip: matlab scripts extract parameters from CESM nc files, 
% regrid the data and save them to .mat file   

% update history:
% v1.0 DL 2019Oct05

% extra notes:
% parameters needed: u,v,w velocity (UVEL,VVEL,WVEL), 
% wind stress (TAUX, TAUY) and temperature (TEMP). 
% =============


%% ====set up environments====
clear all;close all;clc;

  date_str='2019Oct15';
  
%InDir = '/atlantis3/zhangqiuying/regrid_hybrid_v5_rel04_BC5_ne120_t12_pop62_diag.02/ocn_kuroshio/';
InDir1 = '../raw_data/';
WVEL_Str = {'WVEL','time','TLONG','TLAT','z_w_top'}; % WVEL files
TAUX_Str = {'TAUX','time','ULONG','ULAT'};% TAUX files
TAUY_Str = {'TAUY','time','ULONG','ULAT'}; % TAUY files
  
% Infile = [InDir1 'kuroshio_hybrid_v5_rel04_BC5_ne120_t12_pop62_diag.02_WVEL_00870101-00871231.nc'];
% ncdisp(Infile)
%%
[~,ulat,ulon,taux_raw] = LoadCESMKEncFiles3DVarFunc(InDir1,TAUX_Str);
[~,~,~,tauy_raw] = LoadCESMKEncFiles3DVarFunc(InDir1,TAUY_Str);

IndxZ1w = 2;
[time,tlat,tlon,w1_raw,z1_w] = LoadCESMKEncFiles4DVarFunc(InDir1,WVEL_Str,IndxZ1w);

IndxZ2w = 5;
[~,~,~,w2_raw,z2_w] = LoadCESMKEncFiles4DVarFunc(InDir1,WVEL_Str,IndxZ2w);

IndxZ3w = 8;
[~,~,~,w3_raw,z3_w] = LoadCESMKEncFiles4DVarFunc(InDir1,WVEL_Str,IndxZ3w);

IndxZ4w = 15;
[~,~,~,w4_raw,z4_w] = LoadCESMKEncFiles4DVarFunc(InDir1,WVEL_Str,IndxZ4w);

% indxlat = find(tlat(1,:)>33 & tlat(1,:)<35);
% indxlon = find(tlon(:,1)>137 & tlon(:,1)<142);
% w1_raw(indxlon,indxlat,:,:)=nan;
% w2_raw(indxlon,indxlat,:,:)=nan;
% w3_raw(indxlon,indxlat,:,:)=nan;
% w4_raw(indxlon,indxlat,:,:)=nan;

% indx = [1:30 365-29:365];
% w1_av = nanmean(w1_raw(:,:,indx),3);
% w2_av = nanmean(w2_raw(:,:,indx),3);
% w3_av = nanmean(w3_raw(:,:,indx),3);
% w4_av = nanmean(w4_raw(:,:,indx),3);

figure;
  plt_indx = 51;
  indxX = 1 : 5 : 226;
  indxY = 1 : 5 : 401;
quiver(ulon(indxY,indxX),ulat(indxY,indxX), ...
      taux_raw(indxY,indxX,plt_indx),tauy_raw(indxY,indxX,plt_indx), ...
      'b','AutoScale','off');
   hold on;
  quiver(132,44,0.5,0,'b','AutoScale','off');
   text(134,44,'0.5 N/m^2')
   title('wind stress field [N/m^2]') 
   
%% === make pics ===
f1=figure; % check rotation
  set(f1,'units','normalized','position',[0,0,1,1]);
%   50    51   347   364   365 
  plt_indx = 51;
  w1_plt=w1_raw(:,:,plt_indx)./100;
  w1_plt(abs(w1_plt-nanmean(w1_plt(:)))>3*nanstd(w1_plt(:),1))=nan;
  w2_plt=w2_raw(:,:,plt_indx)./100;
  w2_plt(abs(w2_plt-nanmean(w2_plt(:)))>3*nanstd(w2_plt(:),1))=nan; 
  w3_plt=w3_raw(:,:,plt_indx)./100;
  w3_plt(abs(w3_plt-nanmean(w3_plt(:)))>3*nanstd(w3_plt(:),1))=nan;
  w4_plt=w4_raw(:,:,plt_indx)./100;
  w4_plt(abs(w4_plt-nanmean(w4_plt(:)))>3*nanstd(w4_plt(:),1))=nan;   
  
%   w1_plt=w1_av;
%   w2_plt=w2_av;
%   w3_plt=w3_av;
%   w4_plt=w4_av;
  
 subplot(2,2,1) % raw data
  [c,h] = contourf(tlon,tlat,w1_plt);
  set(h,'linestyle','none');colorbar;
  title(['CESM WVel at ', num2str(z1_w./100),' m, indxZ:',num2str(IndxZ1w)]);
 subplot(2,2,2) % data after rotation
  [c,h] = contourf(tlon,tlat,w2_plt);
  set(h,'linestyle','none');colorbar;
  title(['CESM WVel at ', num2str(z2_w./100),' m, indxZ:',num2str(IndxZ2w)]);
 subplot(2,2,3) % data after rotation and interpolation
  [c,h] = contourf(tlon,tlat,w3_plt);
  set(h,'linestyle','none');colorbar;
  title(['CESM WVel at ', num2str(z3_w./100),' m, indxZ:',num2str(IndxZ3w)]);
 subplot(2,2,4) % data after rotation and interpolation
  [c,h] = contourf(tlon,tlat,w4_plt);
  set(h,'linestyle','none');colorbar;
  title(['CESM WVel at ', num2str(z4_w./100),' m, indxZ:',num2str(IndxZ4w)]);
% ==============