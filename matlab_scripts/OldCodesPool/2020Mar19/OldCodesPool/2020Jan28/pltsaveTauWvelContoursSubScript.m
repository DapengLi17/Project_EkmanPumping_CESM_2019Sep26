pic1 = ['../pics/TauWvelContourWinterStorm_KE_',date_str,'.png'];

%% --- Tau, W winter storms ---
f1=figure('Renderer','painters');
  
  % set pic size
  pic1_size=[12 10]; % pic size: unit [inch] 3.155 inch is the width(length in x direction) of JPO 1 column pic, height varies according to your needs
  set(f1,'Units','inches','Position',[5,5,pic1_size]);
  font_size=10;
  ylabel_position=[-0.08, 0.5, 0];
  Xlim     = [130 170];     Xticks=[130:10:170];
  Ylim     = [26.5 45.5];   Yticks=[30:5:45];
  Tau_lim  = [0 0.4]; TauTicks=[0:0.1:0.4];
  W_lim_ce = [-2 2]; W_ce_Ticks =[-2:1:2];
  W_lim_nl = [-2 2]; W_nl_Ticks =[-2:1:2];

% ~~~ generate subplot position ~~~  
  row_num=4;col_num=3;margin_left=0.06;
  margin_right=0.07;margin_top=0.04;margin_botm=0.06;
  pics_dist_x=0.03; pics_dist_y=0.03;
 
  [sbplt_posit]=compute_subplots_position_matrix(row_num,col_num,margin_left, ...
      margin_right,margin_top,margin_botm,pics_dist_x,pics_dist_y);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

 subplot('Position',sbplt_posit(1,:))
  pcolor(Lon,Lat,tau_ns_av);shading interp;caxis(Tau_lim); %hc=colorbar;title(hc,'N/m^2');%polarmap(64,1);
    set(gca,'xlim',Xlim,'XTick',Xticks,'xticklabel',[],'ylim',Ylim, ...
        'ytick',Yticks,'fontsize',font_size);    
    ylabel(['\tau',sprintf('\n'),'Lat'],'Units','Normalized', ...
        'Position',ylabel_position,'fontsize',font_size);
    title(['storm-free days av (N: ',num2str(length(indx_ns) ),' days)']);
    text(132,28,'a','fontsize',font_size,'FontWeight','bold');
    
 subplot('Position',sbplt_posit(2,:))
  pcolor(Lon,Lat,tau_st_av);shading interp;caxis(Tau_lim); %hc=colorbar;title(hc,'N/m^2');%polarmap(64,1);
    set(gca,'xlim',Xlim,'XTick',Xticks,'xticklabel',[],'ylim',Ylim, ...
        'ytick',Yticks,'yticklabel',[],'fontsize',font_size);    
    title(['storm days av (N: ',num2str(length(indx_st) ),' days)']);
    text(132,28,'b','fontsize',font_size,'FontWeight','bold');  
    
 subplot('Position',sbplt_posit(3,:))
  pcolor(Lon,Lat,tau_st_av-tau_ns_av);shading interp;caxis(Tau_lim); %hc=colorbar;title(hc,'N/m^2');%polarmap(64,1);
    set(gca,'xlim',Xlim,'XTick',Xticks,'xticklabel',[],'ylim',Ylim, ...
        'ytick',Yticks,'YTickLabel',[],'fontsize',font_size);    
    title(['storm days av - storm-free days av']);
    text(132,28,'c','fontsize',font_size,'FontWeight','bold');  
    hc=colorbar('FontSize',font_size-2);
    set(hc,'position',[0.96 sbplt_posit(3,2) .01 sbplt_posit(3,4)-0.02], ...
       'ylim',Tau_lim,'YTick',TauTicks);
    title(hc,'[N/m^2]','FontSize',font_size-2);
 
 subplot('Position',sbplt_posit(4,:))
  pcolor(Lon,Lat,w_li_ns_av.*SecPerDay);shading interp;caxis(W_lim_nl); %hc=colorbar;title(hc,'N/m^2');%polarmap(64,1);
    hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',Xlim,'XTick',Xticks,'xticklabel',[],'ylim',Ylim, ...
        'ytick',Yticks,'fontsize',font_size);    
    ylabel(['w_f',sprintf('\n'),'Lat'],'Units','Normalized', ...
        'Position',ylabel_position,'fontsize',font_size);
    text(132,28,'d','fontsize',font_size,'FontWeight','bold');  
  
 subplot('Position',sbplt_posit(5,:))
  pcolor(Lon,Lat,w_li_st_av.*SecPerDay);shading interp;caxis(W_lim_nl); %hc=colorbar;title(hc,'N/m^2');%polarmap(64,1);
    hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',Xlim,'XTick',Xticks,'xticklabel',[],'ylim',Ylim, ...
        'ytick',Yticks,'yticklabel',[],'fontsize',font_size);    
    text(132,28,'e','fontsize',font_size,'FontWeight','bold'); 
    
 subplot('Position',sbplt_posit(6,:))
  pcolor(Lon,Lat,(w_li_st_av-w_li_ns_av).*SecPerDay);shading interp;
    caxis(W_lim_nl);hold on;%hc=colorbar;title(hc,'N/m^2');%polarmap(64,1);
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',Xlim,'XTick',Xticks,'xticklabel',[],'ylim',Ylim, ...
        'ytick',Yticks,'yticklabel',[],'fontsize',font_size);    
    text(132,28,'f','fontsize',font_size,'FontWeight','bold'); 
    hc=colorbar('FontSize',font_size-2);
    set(hc,'position',[0.96 sbplt_posit(6,2) .01 sbplt_posit(6,4)-0.02], ...
       'ylim',W_lim_nl,'YTick',W_nl_Ticks);
    title(hc,'[m/day]','FontSize',font_size-2);
  
 subplot('Position',sbplt_posit(7,:))
  pcolor(Lon,Lat,w_nl_ns_av.*SecPerDay);shading interp;caxis(W_lim_nl); %hc=colorbar;title(hc,'N/m^2');%polarmap(64,1);
    hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',Xlim,'XTick',Xticks,'xticklabel',[],'ylim',Ylim, ...
        'ytick',Yticks,'fontsize',font_size);    
    ylabel(['w_{\xi}',sprintf('\n'),'Lat'],'Units','Normalized', ...
        'Position',ylabel_position,'fontsize',font_size);
    text(132,28,'g','fontsize',font_size,'FontWeight','bold');
    
 subplot('Position',sbplt_posit(8,:))
  pcolor(Lon,Lat,w_nl_st_av.*SecPerDay);shading interp;caxis(W_lim_nl); %hc=colorbar;title(hc,'N/m^2');%polarmap(64,1);
    hold on;
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',Xlim,'XTick',Xticks,'xticklabel',[],'ylim',Ylim, ...
        'ytick',Yticks,'YTickLabel',[],'fontsize',font_size);    
    text(132,28,'h','fontsize',font_size,'FontWeight','bold');
    
 subplot('Position',sbplt_posit(9,:))
  pcolor(Lon,Lat,(w_nl_st_av-w_nl_ns_av).*SecPerDay);
    shading interp;caxis(W_lim_nl);hold on;%hc=colorbar;title(hc,'N/m^2');%polarmap(64,1);
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',Xlim,'XTick',Xticks,'xticklabel',[],'ylim',Ylim, ...
        'ytick',Yticks,'YTickLabel',[],'fontsize',font_size);    
    text(132,28,'i','fontsize',font_size,'FontWeight','bold');
    hc=colorbar('FontSize',font_size-2);
    set(hc,'position',[0.96 sbplt_posit(9,2) .01 sbplt_posit(9,4)-0.02], ...
       'ylim',W_lim_nl,'YTick',W_nl_Ticks);
    title(hc,'[m/day]','FontSize',font_size-2);
    
 subplot('Position',sbplt_posit(10,:))
  pcolor(Lon,Lat,w_ce_ns_av.*SecPerDay);shading interp;caxis(W_lim_ce);
    hold on;%hc=colorbar;title(hc,'N/m^2');%polarmap(64,1);
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',Xlim,'XTick',Xticks,'ylim',Ylim, ...
        'ytick',Yticks,'fontsize',font_size);    
    text(132,28,'i','fontsize',font_size,'FontWeight','bold');
    xlabel('Lon','fontsize',font_size);
    ylabel(['w_{t}',sprintf('\n'),'Lat'],'Units','Normalized', ...
        'Position',ylabel_position,'fontsize',font_size);
    
 subplot('Position',sbplt_posit(11,:))
  pcolor(Lon,Lat,w_ce_st_av.*SecPerDay);shading interp;caxis(W_lim_ce);
    hold on;%hc=colorbar;title(hc,'N/m^2');%polarmap(64,1);
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',Xlim,'XTick',Xticks,'ylim',Ylim, ...
        'ytick',Yticks,'YTickLabel',[],'fontsize',font_size);    
    text(132,28,'j','fontsize',font_size,'FontWeight','bold');
    xlabel('Lon','fontsize',font_size);
   
 subplot('Position',sbplt_posit(12,:))
  pcolor(Lon,Lat,(w_ce_st_av-w_ce_ns_av).*SecPerDay);
    shading interp;caxis(W_lim_ce);hold on;%hc=colorbar;title(hc,'N/m^2');%polarmap(64,1);
  [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
    set(h,'LineWidth',1.5);
    set(gca,'xlim',Xlim,'XTick',Xticks,'ylim',Ylim, ...
        'ytick',Yticks,'YTickLabel',[],'fontsize',font_size);    
    text(132,28,'k','fontsize',font_size,'FontWeight','bold');
    xlabel('Lon','fontsize',font_size);
    hc=colorbar('FontSize',font_size-2);
    set(hc,'position',[0.96 sbplt_posit(12,2) .01 sbplt_posit(12,4)-0.02], ...
       'ylim',W_lim_ce,'YTick',W_ce_Ticks);
    title(hc,'[m/day]','FontSize',font_size-2);
  
    polarmap(64,1);
    
% ~~~ plot Ro ~~~    
%  subplot('Position',sbplt_posit(13,:))
%   pcolor(Lon,Lat,Ro_ns_av);shading interp;
%     caxis(W_lim_nl);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
%      hold on;
%   [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
%     set(h,'LineWidth',1.5);
%     ylabel('Q_{\xi}^{50m}');
%   
%  subplot('Position',sbplt_posit(14,:))
%   pcolor(Lon,Lat,Ro_st_av);shading interp;
%     caxis(Q_lim_nl);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
%      hold on;
%   [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
%     set(h,'LineWidth',1.5);
%     
%  subplot('Position',sbplt_posit(15,:))
%   pcolor(Lon,Lat,Ro_st_av-Ro_ns_av);shading interp;
%     caxis(Q_lim_nl);polarmap(64,1);hc=colorbar;title(hc,'W/m^2');
%      hold on;
%   [c,h]=contour(Lon,Lat,Ro_wt_std,[Ro_std_limit Ro_std_limit],'k');
%     set(h,'LineWidth',1.5);  
% ~~~~~~~~~~~~~~~~~~~~~~~
% ====================


%% === output data ===
print(f1,'-dpng','-r500',pic1)
% ====================
