function [dist indx] = func_NNR_indx_dist(xvec,yvec,x0,y0,NNR_flag)
%
% Description : Function to find the index and distance of Nearest Neighbor
%
% Input : xvec  - column vector of xpoints to search for
%         yvec  -   "              ypoints to search for
%         x0    - scalar x reference point
%         y0    - scalar y reference point
%         NNR_Flag - scalar, float, optional,  1 (default) = Nearest neighbour
%                                              0           = Farthest Neighbour.
%
% Output : indx - index of nearest neighbor in xvec/yvec, with respect to point (x0,y0)
%          distance - distance between (x0,y0) and nearest point in xvec/yvec
%
% Example :  xpts = [1,2,3,4] ;
%            ypts = [11,12,13,14] ;
%            [idx dis] = func_NNR_indx_dist(xpts,ypts,3.5,12) 
%
% Method : dist_square = (x - x0)^2 + (y - y0)^2 
%
% Written By : Jaison Kurian (jaison@atmos.ucla.edu)
% Written On : May/22/2008
%
% June/22/2009 : A flag is added as 5th argument to get the distance and
%                    index for the farthest Neighbor.
% 
%=======================================================================

    % check number of input arguments         


    if (nargin == 4) 
       NNR_flag = 1 ; 
    elseif (nargin == 5)
        if ( NNR_flag ~= 0 && NNR_flag ~= 1 )
            error ('\n\t FATAL (%s) : The 5th argument should be 1 to work with \n\t\t NEAREST Neighbor (default) and 0 to work with FARTHEST Neighbour.\n\n',mfilename)
        end
    else 
        error('func_NNR_indx_dist:argchk', '\n    FATAL (%s) : Too few input arguments. \n',mfilename)
    end

    % find distance and index
    dist1        = sqrt( (xvec - x0).^2 + (yvec - y0).^2 ) ;
    if ( NNR_flag == 1)
       [dist indx]  = min(abs([dist1]')); 
    else 
       [dist indx]  = max(abs([dist1]')); 
    end
  
    clear dist1 

    return
