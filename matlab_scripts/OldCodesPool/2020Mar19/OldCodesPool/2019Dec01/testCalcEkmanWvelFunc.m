clear all;

x = [1 2 3;1 2 3;1 2 3]
% x =
%      1     2     3
%      1     2     3
%      1     2     3

y = [6 6 6;4 4 4;2 2 2]
% y =
%      6     6     6
%      4     4     4
%      2     2     2

taux = [0.3 0.3 0.3;0.2 0.4 0.2;0.1 0.1 0.1]
% taux =
%     0.3000    0.3000    0.3000
%     0.2000    0.4000    0.2000
%     0.1000    0.1000    0.1000

tauy = -[0.1 0.2 0.3;0.4 1.0 0.6;0.7 0.8 0.9]
% tauy =
%    -0.1000   -0.2000   -0.3000
%    -0.4000   -1.0000   -0.6000
%    -0.7000   -0.8000   -0.9000

% test divergence 
div = divergence(x,y,taux,tauy)
% div =
%     0.1500    0.4000    0.1500
%     0.3500    0.1500   -0.0500
%     0.1500   -0.1000    0.1500
    
% left bottom corner
0.3/2
% ans = 0.1500

% center point
0.6/4
% ans = 0.1500

f = [3 3 3;2 2 2;1 1 1]
% f =
%      3     3     3
%      2     2     2
%      1     1     1

kesai = [0.9 0.8 0.7;0.6 0.5 0.4;0.3 0.2 0.1]
% kesai =
%     0.9000    0.8000    0.7000
%     0.6000    0.5000    0.4000
%     0.3000    0.2000    0.1000

% f + kesai
%     3.9000    3.8000    3.7000
%     2.6000    2.5000    2.4000
%     1.3000    1.2000    1.1000
    
rho0 = 2;
[W] = CalcEkmanWvelFunc(rho0,x,y,taux,tauy,f,kesai)

% W =
%    -0.0135    0.0064   -0.0137
%    -0.1231   -0.0235    0.0762
%    -0.0641   -0.0891   -0.0739

% left bottom corner
-(-(-0.8/1.2+0.7/1.3) + (0.2/2.6-0.1/1.3)/2)/rho0 
% ans = 0.0641

% center point
-(-(-0.6/2.4+0.4/2.6)/2 + (0.3/3.8-0.1/1.2)/4)/rho0
% ans = -0.0235 
