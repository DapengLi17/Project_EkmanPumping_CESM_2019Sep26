function [Val,Ind]=min2d(Mat);
%--------------------------------------------------------------------------
% min2d function      2d minimum function. Return the minimum value
%                   and index in a 2d matrix.
% Input  : - Matrix.
% Output : - Minimum value.
%          - [I,J] index of minimum.
% Tested : Matlab 5.3
%     By : Eran O. Ofek                  October 2000
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%--------------------------------------------------------------------------

[V1,I1] = min(Mat);
[V2,I2] = min(V1);

Val = V2;
Ind(1) = I1(I2);   % I
Ind(2) = I2;       % J

