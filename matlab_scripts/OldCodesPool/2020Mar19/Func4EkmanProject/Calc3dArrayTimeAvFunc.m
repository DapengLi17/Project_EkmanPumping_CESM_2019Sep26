
%% ========================  readme  =============================
% 
% DESCRIPTION:
% 
%  A function to compute time average for 3D array (lat,lon,time).
%
% update history:
% v1.0 DL 2019Nov14
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
%   dat_raw   - 3D dat matrix, 1st Dimension: lat, 2nd Dimension: lon, 3rd
%               Dimension: time
%   indx      - time index, mostly commonly for month or season, it can be 
%               computed from FindMonthlyTimeIndxFunc.m
%   av_method - string for average method, either 'nanmean' or 'bootstrap'
% 
% OUTPUT:
%   dat_av    - 2D dat matrix, averaged over time based on the index 
%
% EXTRA NOTES:
% N/A 
% 
% REFERENCE:
% N/A
% ====================================================================

function [dat_av] = Calc3dArrayTimeAvFunc(dat_raw,indx,av_method)

%% === data analysis ===
m   = size(dat_raw,1);
n   = size(dat_raw,2);
dat = dat_raw(:,:,indx);

if strcmp(av_method,'bootstrap');
 for i = 1 : m
  for j = 1 : n
   dat1d = dat(i,j,:);
   [~,x_mn,~] = bootstrap5(dat1d(:));
   dat_av(i,j) = x_mn;
   clear dat1d x_mn
  end
 end
elseif strcmp(av_method,'nanmean');
 dat_av = nanmean(dat,3);   
else
 disp('av methods must be bootstrap or nanmean !')    
end
% ======================

end
