function [ubar,vbar]=func_roms_geovel(f,zeta,pm,pn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2000 IRD                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                 %
%                                                                 %
%   compute the geostrophy...                                     % 
% copy of surf_geostr_vel of IRD Roms_Tools                       %
%                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=9.81;
[Lp,Mp]=size(zeta);
L=Lp-1;
M=Mp-1;
Lm=L-1;
Mm=M-1;
ubar=zeros(Lp,M);
vbar=zeros(L,Mp);
zetap=zeros(L,M);
zetap=(zeta(1:L,1:M)+zeta(2:Lp,1:M)+...
       zeta(2:Lp,2:Mp)+zeta(1:L,2:Mp))/4;
vbar(1:L,2:M)=g.*(pn(1:L,2:M)+pn(2:Lp,2:M)).*...
              (zetap(1:L,2:M)-zetap(1:L,1:Mm))./...
              (f(1:L,2:M)+f(2:Lp,2:M));
vbar(:,1)=vbar(:,2);
vbar(:,Mp)=vbar(:,M);
ubar(2:L,1:M)=-g.*(pm(2:L,1:M)+pm(2:L,2:Mp)).*...
              (zetap(2:L,1:M)-zetap(1:Lm,1:M))./...
              (f(2:L,1:M)+f(2:L,2:Mp));
ubar(1,:)=ubar(2,:);
ubar(Lp,:)=ubar(L,:);
return

