x = [1 2;1 2];
% x =
%      1     2
%      1     2
y = [2 2;1 1];
% y =
%      2     2
%      1     1

taux = [0.1 0.2;0.3 0.4];
% taux =
%     0.1000    0.2000
%     0.3000    0.4000
tauy = [1.1 1.2;1.3 1.4];
% tauy =
%     1.1000    1.2000
%     1.3000    1.4000

W1 = curl(x,y,taux,tauy)
% W1 =
%     0.1500    0.1500
%     0.1500    0.1500
W2 = -divergence(x,y,-tauy,taux)
% W2 =
%     0.3000    0.3000
%     0.3000    0.3000
