function z = func_zlev(h,theta_s,theta_b,Tcline,N,type);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                
% Original name of the function was zlev.
%                                                           %
% function z = zlev(h,theta_s,theta_b,Tcline,N,type);                       %             %
%                                                                           %
% this function compute the depth of rho or w points                        %
%                                                                           %
% On Input:                                                                 %
%                                                                           %
%    type    'r': rho point 'w': w point                                    %
%                                                                           %
% On Output:                                                                %
%                                                                           %
%    z       Depths (m) of RHO- or W-points (3D matrix).                    %
%                                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<6, type='r'; end;

hmin=min(min(h));
hc=min(hmin,Tcline);
[M,L]=size(h);
ds=1/N;
zeta(1:M,1:L)=0;

%----------------------------------------------------------------------------
% Define S-Curves at vertical RHO- or W-points (-1 < sc < 0).
%----------------------------------------------------------------------------

if type=='w',
  sc=-1+(0:N)*ds;
  N=N+1;
else
  sc=-1+((1:N)-0.5)*ds;
end;

%disp(num2str(size(theta_b)))
%disp(num2str(size(theta_s)))
%disp(num2str(size(sc)))

Ptheta=sinh(theta_s.*sc)./sinh(theta_s);
Rtheta=tanh(theta_s.*(sc+0.5))./(2*tanh(0.5*theta_s))-0.5;

%disp(num2str(size(Ptheta)))
%disp(num2str(size(Rtheta)))

Cs=(1-theta_b).*Ptheta+theta_b.*Rtheta;

z=zeros(N,M,L);
for k=1:N
  z(k,:,:)=zeta.*(1+sc(k))+hc.*sc(k)+ ...
           (h-hc).*Cs(k);
end
return
