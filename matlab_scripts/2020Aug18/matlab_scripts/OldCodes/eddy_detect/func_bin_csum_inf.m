function [csum bins bstep count] = func_bin_csum_inf (field, Nbins, bmin, bmax)
%
% Description : Matlab function to bin and find cumulative sum in the reverse
%               order, to make Y-log plots of eddy counts etc.
%
%               say the counts are   1   2   3   4  5
%               then, csum will be   15  14  12  9  5 
%
%                
% Reference : Script from Francois Colas and func_pdf.m
%
% Inputs  : 1  field (float, Matrix, 2D-4D) Variable to find PDF
%           2  Nbins (float, scalar, optional) number of bins 
%           3  bmin  (float, scalar, optional) minimum value for binning
%           4  bmax  (float, scalar, optional) maximum value for binning
%
% Outputs : 1  csum (float, 1D-vector) cumulative sum of bins
%           2  bins  (float, 1D-vector) Bins
%           3  bstep (float, scalar)    spacing between bins
%           4  count (flaot, 1D-vector) original (not cumulative) sum
%
% Usage  : [count bins bstep] = func_pdf(field, 100)
%
%          histogram     : bar(bins,count) 
%                          bar(bins,count/bstep)  ! to get the area=1 
%          line/gaussian : plot(bins, count)
%
%          log-y         : figure ; semilogy(bins, count) ; grid on
 
%
% Jaison Kurian, Francois Colas
% IGPP UCLA 
%           Sep/10/2009
%           Sep/24/2009, bins are not modified now. It make sense to plot
%                           all values of say....10 km and above against
%                           10 km itself, instead of mid-value of that bin-class
%                           say...15 or 20 km.
%           Oct/04/2010, Now, the end of Bins is +inf, manually done within this 
%                           function.  
%----------------------------------------------------------------------------
 
% Initialize

     if ( nargin < 1 ) 
         error ('\n\t FATAL (%s) : Need a variable to find PDF.\n', mfilename)
     elseif ( nargin == 1)
         Nbins = 20 ;
     elseif ( nargin > 1 )
         Nbins = abs(Nbins) ;
     end          

% make arrays to 1D-vector


     N        = numel(field) ;           % total number of points, including NaN
     field    = reshape(field,1,N) ;     % reshape to 1D-vector 

% set bins
     
     if ( nargin < 3 )
         bmax = max(abs(field)) ;
         bmin = -1 * bmax ;
     elseif ( nargin == 3 )
         bmax = max(abs(field)) ;
     end 

     bstep = (bmax-bmin)/(Nbins-1) ;
     bins  = [bmin:bstep:bmax] ;

% find histograms

     bins  = [bins +inf] ;
     count = histc(field, bins) ;
% find cumulative sum in the required direction

     csum = zeros(size(count)) ;
     for i = 1:length(csum) 
        csum(i) = sum(count(i:end)) ;
     end     

% DONE

   
     return
