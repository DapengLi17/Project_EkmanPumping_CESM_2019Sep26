%% ===readme===

% descrip: matlab scripts plot vertical velocities from 
% 1) CESM model output
% 2) Stern nonlinear Ekman vertical pumping velocity
% 3) linear Ekman vertical pumping velocity

% update history:
% v1.0 DL 2019Oct07

% extra notes:
% =============


%% ====set up environments====
clear all;close all;clc;

  date_str='2019Oct16';
  
% InFile1 = '../data_after_manipulation/CESMtauuvwtempKE_2019Oct04.mat';
InFile1 = '../data_after_manipulation/CESMtauuvwtempKE_2019Oct15.mat';
% Pic1 = ['../pics/WindStress_',date_str,'.png'];
% Pic2 = ['../pics/SfcCurrent_',date_str,'.png'];
% Pic3 = ['../pics/Vorticity_',date_str,'.png'];
% Pic4 = ['../pics/CESM_Wvel_',date_str,'.png'];
% Pic5 = ['../pics/Stern_Wvel_',date_str,'.png'];

Pic1 = ['../pics/CESM_TauUVvelKesaiRoNum_',date_str,'.png'];
Pic2 = ['../pics/CESM_WvelEkmanvel_',date_str,'.png'];

% av_method='bootstrap';
% addpath('/home/dapengli/MatlabCodes4DatAnalysis_DL/added_Matlab_toolbox_for_Research_30Aug18/smooth2a')
addpath(genpath('/home/dapengli/MatlabCodes4DatAnalysis_DL/'))
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
w = KE.w;
taux = KE.taux;
tauy = KE.tauy;
f = repmat(sw_f(KE.Lat_vec),1,length(KE.Lon_vec));

for i = 1 : size(u,3) % loop through time
    kesai(:,:,i) = CalcRelVorticityFunc(x,y,u(:,:,i),v(:,:,i));
end

for i = 1 : size(kesai,3) % loop through time
 [W_NL(:,:,i)] = CalcEkmanWvelFunc(rho0,x,y, ...
     taux(:,:,i),tauy(:,:,i),f,kesai(:,:,i));
end

kesai0 = zeros(size(x));
for i = 1 : size(kesai,3) % loop through time
 [W_LI(:,:,i)] = CalcEkmanWvelFunc(rho0,x,y, ...
     taux(:,:,i),tauy(:,:,i),f,kesai0);
end

% TotVor = f + kesai; % total vorticity
for i = 1 : size(kesai,3)
    Ro(:,:,i) = kesai(:,:,i)./f; % Rossby number 
end

W_NL_cln=W_NL;
W_LI_cln=W_LI;
  
% W_NL_cln(abs(TotVor)<5e-5)=nan;
% W_LI_cln(abs(TotVor)<5e-5)=nan;

% W_NL_cln(abs(Ro)>0.5)=nan;
% W_LI_cln(abs(Ro)>0.5)=nan;

% --- check Ekman W vel ---
% [N_W_NL,C_W_NL] = hist(W_NL_cln(:),100); 
% [N_W_LI,C_W_LI] = hist(W_LI_cln(:),100); 
% figure;
% h_Np=bar(C_W_NL,N_W_NL,'b','FaceAlpha',0.2, ...
%       'EdgeAlpha',.8,'linewidth',1.2);hold on;
% h_Sp=bar(C_W_LI,N_W_LI,'r','FaceAlpha',0.2, ...
%       'EdgeAlpha',.8,'linewidth',1.2);grid on;

% indxX = [1:3];
% indxY = [length(y_vec)-2:length(y_vec)]
% f(indxY,indxX)
% x(indxY,indxX)
% y(indxY,indxX)
% taux(indxY,indxX,1)
% tauy(indxY,indxX,1)
% 
% -tauy(indxY,indxX,1)./f(indxY,indxX)
% taux(indxY,indxX,1)./f(indxY,indxX)
% 
% W_LI(indxY,indxX,1) % compare the center point of the output matrix
% ((1.950-1.909)*10^3/1.98/1e4 + (-81.23+151.23)/2.2/1e4)/(-rho0)
% ================


%% === make pics ===
f1=figure; % surface wind

  plt_indx = 51;
  indxX = 1 : 5 : length(KE.Lon_vec);
  indxY = 1 : 5 : length(KE.Lat_vec);
  set(f1,'units','normalized','position',[0 0 1 1])
  
 subplot(2,2,1)
  quiver(KE.Lon(indxY,indxX),KE.Lat(indxY,indxX), ...
      taux(indxY,indxX,plt_indx),tauy(indxY,indxX,plt_indx), ...
      'b','AutoScale','off');
   hold on;
  quiver(132,44,0.5,0,'b','AutoScale','off');
   text(134,44,'0.5 N/m^2')
   ylabel('Lat [\circ]')
   title('wind stress field [N/m^2] wind spd ~ 25 m/s') 
   
  subplot(2,2,2) 
   quiver(KE.Lon(indxY,indxX),KE.Lat(indxY,indxX), ...
      u(indxY,indxX,plt_indx),v(indxY,indxX,plt_indx), ...
      'b','AutoScale','off');hold on;
   quiver(132,44,1,0,'b','AutoScale','off');
   set(gca,'xlim',[130 170],'ylim',[26 46]);
   text(134,44,'1 m/s');
   title(['CESM horizontal vel [m/s] dpth:',num2str(KE.dpth_u),'m']);
   
 subplot(2,2,3)
  [c,h] = contourf(KE.Lon,KE.Lat,kesai(:,:,plt_indx));
  set(h,'LineStyle','none');colorbar;
  caxis([-6 6]*1e-5);polarmap(64,1)
  set(gca,'YLim', [27 45],'ytick',[27:3:45]);
  ylabel('Lat [\circ]');xlabel('Lon [\circ]');
  title(['CESM relative vorticity [s^{-1}] dpth:',num2str(KE.dpth_u),'m']);
  
 subplot(2,2,4)
  [c,h] = contourf(KE.Lon,KE.Lat,kesai(:,:,plt_indx)./f);colorbar
  set(h,'LineStyle','none');colorbar;
  caxis([-1 1]);polarmap(64,1)
%   hold on;
%   [c,h] = contour(KE.Lon,KE.Lat,abs(kesai(:,:,plt_indx)./f),[0.5 0.5],'r');
  title(['CESM Rossby number dpth:',num2str(KE.dpth_u),'m']);
  xlabel('Lon [\circ]');

%%  
f2=figure; 
  set(f2,'units','normalized','position',[0 0 1 1])  
 subplot(2,2,1) % CESM W vel with u v quiver
  w_cesm_plt=w(:,:,plt_indx);
  w_cesm_plt(abs(w_cesm_plt-nanmean(w_cesm_plt(:)))>3*nanstd(w_cesm_plt(:),1))=nan;
  Nr=10;Nc=10;
  w_cesmsmo_plt = smooth2a(w_cesm_plt,Nr,Nc);
  
  [c,h] = contourf(KE.Lon,KE.Lat,w_cesm_plt);
%   pcolor(KE.Lon,KE.Lat,w_cesmsmo_plt);shading interp;
  set(h,'linestyle','none');colorbar;
  caxis([-1.5 1.5].*1e-4);polarmap(64,1);
  set(gca,'YLim', [27 45],'ytick',[27:3:45]);
  ylabel('Lat [\circ]')
  title(['CESM raw Wvel [m/s], dpth:',num2str(KE.dpth_w),' m'])
%   hold on
%   quiver(KE.Lon(indxY,indxX),KE.Lat(indxY,indxX), ...
%       u(indxY,indxX,plt_indx),v(indxY,indxX,plt_indx), ...
%       'b','AutoScale','off');
%    quiver(132,44,1,0,'b','AutoScale','off');
%    text(134,44,'1 m/s')
%    title('CESM W (color) and U (quiver) [m/s]')   

 subplot(2,2,2) % linear Ekman vertical pumping velocity 
  [c,h] = contourf(KE.Lon,KE.Lat,w_cesmsmo_plt);
%   [c,h] = contourf(KE.Lon,KE.Lat,W_LI(:,:,plt_indx));
%   [c,h] = contourf(KE.Lon,KE.Lat,W_LI);
  caxis([-1.5 1.5]*1e-4);
  set(h,'linestyle','none');colorbar;
  set(gca,'YLim', [27 45],'ytick',[27:3:45]);
%   ylabel('Lat [\circ]')  
  title(['CESM 1 \circ smo Wvel [m/s], dpth:',num2str(KE.dpth_w),' m'])
%   title('raw Ekman (no \xi) vel [m/s]')   
   
 subplot(2,2,3) % nonlinear Ekman vertical pumping velocity (cut vel > 5e-5 m/s)
  W_NL_plt=W_NL(:,:,plt_indx);
  W_NL_plt = smooth2a(W_NL_plt,Nr,Nc);
%   pcolor(KE.Lon,KE.Lat,W_NL_plt);shading interp;
  [c,h] = contourf(KE.Lon,KE.Lat,W_NL_plt);
%   [c,h] = contourf(KE.Lon,KE.Lat,mean(W_NL,3));
  caxis([-1.5 1.5]*1e-4);
  set(h,'linestyle','none');colorbar;
  set(gca,'YLim', [27 45],'ytick',[27:3:45]);
  ylabel('Lat [\circ]');
  title(['nonlinear Ekman Wvel (Stern 1965) [m/s],dpth:',num2str(KE.dpth_w),' m'])
  
 subplot(2,2,4) % linear Ekman vertical pumping velocity 
  W_LI_plt=W_LI_cln(:,:,plt_indx);
  W_LI_plt = smooth2a(W_LI_plt,Nr,Nc);
  [c,h] = contourf(KE.Lon,KE.Lat,W_LI_plt);
  caxis([-1.5 1.5]*1e-4);set(h,'linestyle','none');colorbar;
  set(gca,'YLim', [27 45],'ytick',[27:3:45]);
  xlabel('Lon [\circ]')
  title(['linear Ekman Wvel [m/s],dpth:',num2str(KE.dpth_w),' m'])
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
print(f2,'-dpng',Pic2)
% print(f3,'-dpng',Pic3)
% print(f4,'-dpng',Pic4)
% print(f5,'-dpng',Pic5)
% ==================