function [xmin,xmean,xmax] = bootstrap5(x,bootnum,ci);

% THIS WORKS ON MATLAB 5.2.1
% This function calculates the mean (xmean) and bootstrap
% confidence points xmin and xmax for an input
% vector of measured values (x).  
% One chooses bootnum as the number of bootstrap samples 
% (500-1000 is reasonable). The optional confidence interval 
% parameter, ci, should be input as fraction, i.e. for a 95% 
% confidence interval enter 0.95 (the default value is 0.95, and 
% will be used if ci is not given as an input).
%
% Original coding: Blair Greenan
% Changes: added option for confidence interval as an input argument (Ayal Anis)

% Get rid of NaNs
index_finite = find(isfinite(x));
x            = x(index_finite);

if nargin == 1 
   bootnum=500;
   ci = 0.95;
elseif nargin == 2 
   ci = 0.95;
end

if (ci > 1) 
  disp(' ')
  disp(' NOTE: confidence interval should be < 1')
  disp(' ')
  return
end

N = length(x);
rand_matrix = 1 + round((N-1)*rand(bootnum,N));
x_resampled = zeros(bootnum,N);
x_resampled = x(rand_matrix);

bstrap_samples = mean(x_resampled,2);
bootlow=floor(bootnum*(1-ci)/2);
boothigh=ceil(bootnum*(ci+(1-ci)/2));

sort_smpl=sort(bstrap_samples);
xmin=sort_smpl(bootlow);
xmean=mean(bstrap_samples);
xmax=sort_smpl(boothigh);
