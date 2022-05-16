%% === readme ===

% c = 'r' or 'b'
% x_min,x_av,x_max,y are column vec
% h_fil: handle of fill
% h_plt: handle of plot
% ===============

function [h_fil h_plt]=plt_shade_line(x_min,x_av,x_max,y,c)
 
 h_fil=fill([x_min(1);x_max;flipud(x_min)],[y(1);y;flipud(y)],c);hold on;
  set(h_fil,'FaceColor',c,'EdgeAlpha',0.0,'FaceAlpha',0.2);
 h_plt=plot(x_av,y,c,'linewidth',1.5);

end

