%% ===readme===

% descrip: matlab scripts plot f and relative vorticity  

% update history:
% v1.0 DL 2019Oct07

% extra notes:
% =============


%% ====set up environments====
clear all;close all;clc;

  date_str='2019Oct07';
  
InFile1 = '../data_after_manipulation/CESMtauuvwtempKE_2019Oct04.mat';
Pic1 = ['../pics/fRelVorticity_',date_str,'.png'];

% av_method='bootstrap';
addpath('seawater_v3.3.1') 
rho0 = 1020;
%=============================


%% === load data ===
KE = load(InFile1); % KE : Kuroshio Extension
%===================


%% === data analysis ===
% convert lat and lon to x and y
% Lat_vec=KE.Lat_vec,Lon_vec=KE.Lon_vec
[x_vec,y_vec,x,y] = LatLon2XYFunc(KE.Lon_vec,KE.Lat_vec);

% compute relative vorticity
u = KE.u;
v = KE.v;
for i = 1 : size(u,3) % loop through time
    kesai(:,:,i) = CalcRelVorticityFunc(x,y,u(:,:,i),v(:,:,i));
end
% ================


%% === make pics ===
f1=figure;

   set(f1,'units','normalized','position',[0 0 1 1])
   plt_indx = 800;
 
  yyaxis left
  [c,h] = contourf(KE.Lon,KE.Lat,kesai(:,:,plt_indx));
  caxis([-8e-5 8e-5]);polarmap(64,2);colorbar('northoutside');
  set(gca,'YLim', [27 45],'ytick',[27:3:45]);
  ylabel('Lat [\circ]')
  yyaxis right
  set(gca,'YLim', [27 45],'ytick',[27:3:45],'YTickLabel',num2str(sw_f(27:3:45),'%7.1e\n'));
  ylabel('f [1/s]')
% ======================

% --- average the whole field through time would smear out eddy 
% (positive + negative = 0), so I do not av through time --- 
% jultime_vec = KE.jultime_vec(1:end-1,:);
% 
% for i = 1 : 12
%   indx = find (jultime_vec(:,2) == i);
%   kesai_dummy = kesai(:,:,indx);
%   kesai_mn(:,:,i) = mean(kesai_dummy,3);
%   clear indx kesai_dummy
% end
% 
% f1=figure;
%  set(f1,'units','normalized','position',[0 0 1 1])
%  for i = 1 : 12
%    subplot(4,3,i)
%      yyaxis left
%   [c,h] = contourf(KE.Lon,KE.Lat,kesai_mn(:,:,i));
%   caxis([-8e-5 8e-5])
%  %colorbar('northoutside');
%   set(gca,'YLim', [27 45],'ytick',[27:3:45]);
%   ylabel('Lat [\circ]')
%   yyaxis right
%   set(gca,'YLim', [27 45],'ytick',[27:3:45],'YTickLabel',num2str(sw_f(27:3:45),'%7.1e\n'));
%   ylabel('f [1/s]')
%  end
% ======================


%% === save data ===
print(f1,'-dpng',Pic1)
% ==================