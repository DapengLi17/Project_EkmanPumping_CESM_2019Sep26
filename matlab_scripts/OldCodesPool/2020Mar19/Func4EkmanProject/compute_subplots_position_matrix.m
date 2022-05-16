%% ===readme===

% descrip: matlab function compute subplots position vectors 

% update history:
% v1.0 DL 2017May07
% v1.1 DL 2018Sep23

% extra notes:
% you need to run this function several times to find the best subplot 
% position for your picture, it depends on font size, page size, etc. 

% input: 
% Nrow: number of rows
% Ncol: number of columns
% LM: left margin [unit is normalized from 0-1]
% RM: right margin [unit is normalized from 0-1]
% TM: top margin [unit is normalized from 0-1]
% BM: bottom margin [unit is normalized from 0-1]
% dist_x: picture distane in x direction (if only 1 column, pics_dist_x=0)
% dist_y: picture distane in y direction (if only 1 row, pics_dist_y=0)

% output
% subplt_posit: subplot position matrix [unit is normalized from 0-1]
% subplt_posit(1,:) is top left picture position vector 
% subplt_posit(2,:) is top right picture position vector (if any)
% subplt_posit(3,:) is bottom left picture position vector (if any)
% subplt_posit(4,:) is bottom right picture position vector (if any)
% =============

function [subplt_posit]=compute_subplots_position_matrix(Nrow,Ncol,LM,RM,TM,BM,dist_x,dist_y)

%% ====data analyzation=======
subplt_posit=nan.*ones(Nrow.*Ncol,4);

% pic length in x axis [unit is normalized from 0-1]
len_x=(1-LM-RM-(Ncol-1).*dist_x)./Ncol;  
% pic length in y axis [unit is normalized from 0-1]
len_y=(1-TM-BM-(Nrow-1).*dist_y)./Nrow;  

counter = 1 ;
 for i = 1 : Nrow % loop through number of rows
    
    for j = 1 : Ncol % loop through number of columns

        left=LM+(j-1)*len_x+(j-1)*dist_x;
        bottom=1-TM-i*len_y-(i-1)*dist_y;
        subplt_posit(counter,:)=[left bottom len_x len_y];
        
        clear left bottom 
        counter=counter+1;
            
    end

 end
% ============================

end
