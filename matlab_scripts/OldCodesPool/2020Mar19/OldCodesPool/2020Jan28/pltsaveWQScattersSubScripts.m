%% --- W,Q scatters ---
f3 = figure('Renderer','painters'); % scatter plot
 
   % set pic size
  pic3_size=[10 7]; % pic size: unit [inch] 3.155 inch is the width(length in x direction) of JPO 1 column pic, height varies according to your needs
  set(f3,'Units','inches','Position',[5,5,pic3_size]);
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
   plot(w_nl_df_Ro,w_li_df_Ro,'b*');hold on;
   plot(W_lim,bW_nl2li(1)+bW_nl2li(2).*W_lim,'r','linewidth',1.5);grid on;
    set(gca,'xlim',[-4e-5 8e-5],'ylim',[-4e-5 6e-5],'fontsize',font_size);
    xlabel('W_{\xi} [m/s]','fontsize',font_size);
    ylabel('W_f [m/s]','fontsize',font_size);
    title(['y=',num2str(bW_nl2li(1),'%.2e'),'+',num2str(bW_nl2li(2),'%.2f'), ...
        '*x, r = ', num2str(rW_nl2li,'%.2f')])
%   [x_min,x_mn,x_max]=bootstrap5(abs(w_nl_Ro./w_li_Ro))
  
  subplot('position',sbplt_posit(2,:))
   plot(w_nl_df_Ro,w_ce_df_Ro,'b*');hold on;
   plot(W_lim,bW_nl2ce(1)+bW_nl2ce(2).*W_lim,'r','linewidth',1.5);grid on;
    set(gca,'xlim',[-4e-5 8e-5],'ylim',[-4e-4 4e-4],'ytick',[-4:2:4].*1e-4, ...
        'fontsize',font_size);
    xlabel('W_{\xi} [m/s]','fontsize',font_size);
    ylabel('W_{cesm} [m/s]','fontsize',font_size);
    title(['y=',num2str(bW_nl2ce(1),'%.2e'),'+',num2str(bW_nl2ce(2),'%.2f'), ...
        '*x, r = ',num2str(rW_nl2ce,'%.2f')])

  subplot('position',sbplt_posit(3,:))
   plot(Q_nl_df_Ro,Q_li_df_Ro,'b*');hold on;
   plot(Q_lim,bQ_nl2li(1)+bQ_nl2li(2).*Q_lim,'r','linewidth',1.5);grid on;
    set(gca,'xlim',[-200 400],'ylim',[-150 150],'fontsize',font_size)
    xlabel('Q_{\xi} [W/m^2]','fontsize',font_size);
    ylabel('Q_f [W/m^2]','fontsize',font_size);
    title(['y=',num2str(bQ_nl2li(1),'%.2f'),'+',num2str(bQ_nl2li(2),'%.2f'), ...
        '*x, r = ',num2str(rQ_nl2li,'%.2f')])
    
  subplot('position',sbplt_posit(4,:))
   plot(Q_nl_df_Ro_Dspk,Q_ce_df_Ro_Dspk,'b*');hold on;
   plot(Q_lim,bQ_nl2ce(1)+bQ_nl2ce(2).*Q_lim,'r','linewidth',1.5);grid on;
    set(gca,'xlim',[-200 400],'ylim',[-500 1500],'fontsize',font_size)
    xlabel('Q_{\xi} [W/m^2]','fontsize',font_size);
    ylabel('Q_{cesm} [W/m^2]','fontsize',font_size);
    title(['y=',num2str(bQ_nl2ce(1),'%.2f'),'+',num2str(bQ_nl2ce(2),'%.2f'), ...
        '*x, r = ',num2str(rQ_nl2ce,'%.2f')])
% --------------------------
