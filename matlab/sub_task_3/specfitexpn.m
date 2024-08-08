function [h,y,resl,rsq] = specfitexpn(param,xdat,ydat)
%param is the list of initial guess parameters
%xdat is the x axis of data to be fit
%ydat is the data to be fit
%
%h is the list of output parameters
%fity is a the fitted curve on xdat
%res is the relative residual between the fit and the data at each point

x0=param;
if length(param)==3
xL=[0 param(2)*0.01 0];
xH=[param(1)*2+0.00001 param(2)*1.05+0.00001 param(3)*3+0.00001];
end

if length(param)==5
xL=[0 param(2)*0.01 0 param(4)*0.0001 param(5)*0.0000001];
xH=[param(1)*2+0.00001 param(2)*1.05+0.00001 param(3)*3+0.00001 param(4)*10+0.00001 param(5)*1.01+0.00001];
end

if length(param)==7
xL=[-param(1)*2-0.0001 0.001 -param(3)*2-0.0001  0.001 -param(5)*2-0.0001  0.001 -0.1];
xH=[param(1)*2+0.00001 param(2)*1.05+0.00001 param(3)*3+0.00001 param(4)*10+0.00001 param(5)*2+0.00001 param(6)*1.01+0.00001 0.1];
end

options = optimset('TolFun',10^-15,'TolX',10^-25,'MaxFunEvals',2000,'MaxIter',10000,'FinDiffType','central');
[h]=lsqcurvefit('expnfitt',x0,xdat,ydat',xL,xH,options);

%Generate curve of final fit
t=xdat;

[y]=expnfitt(h,t);

for z=1:length(t)
resl(z,1)=(y(z)-ydat(z))/mean(ydat);
resi(z,1)=(y(z)-mean(ydat))^2;
toti(z,1)=(ydat(z)-mean(ydat))^2;
end

SSres=sum(resi);
SStot=sum(toti);
rsq=1-(SSres/SStot);

return






