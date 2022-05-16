function [Val,Ind]=max2d(Mat);
%--------------------------------------------------------------------------
% max2d function      2d maximum function. Return the maximum value
%                   and index in a 2d matrix.
% Input  : - Matrix.
% Output : - Maximum value.
%          - [row,column] index of maximum.
% Based on Eran O. Ofek's min2d function
% Tested : Matlab 5.3
%     By : Eran O. Ofek                  October 2000
%    URL : http://wise-obs.tau.ac.il/~eran/matlab.html
%
% Written by : Jaison Kurian
% Written On : June/20/2008
%--------------------------------------------------------------------------

[V1,I1] = max(Mat);
[V2,I2] = max(V1);

Val = V2;
Ind(1) = I1(I2);   % row
Ind(2) = I2;       % column

