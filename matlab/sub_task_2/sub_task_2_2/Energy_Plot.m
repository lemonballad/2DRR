clear all
load EQ_HP2;

h_static_range=2:47;
h_dynamic_range=[1 48:49];
p_static_range=[1:1537];
eps0=8.85*10^-12;
Ang_per_m=10^10;C_per_e=1.6*10^-19;
wvn_per_J=5.034*10^22;wvn_mol_per_kJ=83.593;
prefactor=Ang_per_m*C_per_e^2*wvn_per_J/4/pi/eps0;

r0=10*[hx(1) hy(1) hz(1)];
dr=0.01;
r=-1:dr:1;
lr=length(r);

px=10*px-r0(1);
py=10*py-r0(2);
pz=10*pz-r0(3);
hx=10*hx-r0(1);
hy=10*hy-r0(2);
hz=10*hz-r0(3);
r0=[sum(hx(2:5)) sum(hy(2:5)) sum(hz(2:5))]/4;

v1=[hx(2)-hx(3) hy(2)-hy(3) hz(2)-hz(3)];
v2=[hx(2)-hx(4) hy(2)-hy(4) hz(2)-hz(4)];
n=cross(v1,v2);n=n/norm(n);

EC(1:lr)=0;
LJ(1:lr)=0;

[HX,PX]=ndgrid(hx(h_static_range),px);
[HY,PY]=ndgrid(hy(h_static_range),py);
[HZ,PZ]=ndgrid(hz(h_static_range),pz);
[HQ,PQ]=ndgrid(hq(h_static_range),pq);
[H6,P6]=ndgrid(h6(h_static_range),p6);
[H12,P12]=ndgrid(h12(h_static_range),p12);

Q=HQ.*PQ;
C6=sqrt(H6.*P6);
C12=sqrt(H12.*P12);
XHP=HX-PX;
YHP=HY-PY;
ZHP=HZ-PZ;
R=sqrt((XHP).^2+(YHP).^2+(ZHP).^2);

% EC(1:lr)=sum(sum(Q./R));
% LJ(1:lr)=sum(sum(C12./(R.^12)-C6./(R.^6)));
% [EC(1) LJ(1)]
% EC(1:lr)=prefactor*EC(1:lr);
% LJ(1:lr)=wvn_mol_per_kJ*LJ(1:lr);
% [EC(1) LJ(1)]
% 
[HQ,PQ]=ndgrid(hq(h_dynamic_range),pq(p_static_range));
[H6,P6]=ndgrid(h6(h_dynamic_range),p6(p_static_range));
[H12,P12]=ndgrid(h12(h_dynamic_range),p12(p_static_range));
Q=HQ.*PQ;
C6=sqrt(H6.*P6);
C12=sqrt(H12.*P12);

for ii=1:lr
    [HX,PX]=ndgrid(hx(h_dynamic_range)+r(ii)*n(1),px(p_static_range));
    [HY,PY]=ndgrid(hy(h_dynamic_range)+r(ii)*n(2),py(p_static_range));
    [HZ,PZ]=ndgrid(hz(h_dynamic_range)+r(ii)*n(3),pz(p_static_range));
    
    R=sqrt((HX-PX).^2+(HY-PY).^2+(HZ-PZ).^2);
    
    EC(ii)=EC(ii)+prefactor*sum(sum(Q./R));
    LJ(ii)=LJ(ii)+wvn_mol_per_kJ*sum(sum(C12./R.^12-C6./R.^6));
end

[HQ,PQ]=ndgrid(hq(h_dynamic_range),hq(h_static_range));
[H6,P6]=ndgrid(h6(h_dynamic_range),h6(h_static_range));
[H12,P12]=ndgrid(h12(h_dynamic_range),h12(h_static_range));
Q=HQ.*PQ;
C6=sqrt(H6.*P6);
C12=sqrt(H12.*P12);

for ii=1:lr
    [HX,PX]=ndgrid(hx(h_dynamic_range)+r(ii)*n(1),hx(h_static_range)+r(ii)*n(1));
    [HY,PZ]=ndgrid(hy(h_dynamic_range)+r(ii)*n(2),hy(h_static_range)+r(ii)*n(2));
    [HZ,PY]=ndgrid(hz(h_dynamic_range)+r(ii)*n(3),hz(h_static_range)+r(ii)*n(3));
    
    R=sqrt((HX-PX).^2+(HY-PY).^2+(HZ-PZ).^2);
    
    EC(ii)=EC(ii)+prefactor*sum(sum(Q./R));
    LJ(ii)=LJ(ii)+wvn_mol_per_kJ*sum(sum(C12./R.^12-C6./R.^6));
end

ESS=EC+LJ;

labfont=16;lwidax=2;
figure;plot(r-norm(r0),LJ,'b-','linewidth',4);
hold on
plot([-0.35 -0.35],[min(LJ) max(LJ)],'k:','linewidth',4);
% plot([-2.19 -2.19],[min(LJ) max(LJ)],'r:','linewidth',4);
% plot([-2.1 -2.1],[min(LJ) max(LJ)],'k:','linewidth',4);
plot([-0.24 -0.24],[min(LJ) max(LJ)],'r:','linewidth',4);
axis square;xlabel('dr (Angs)');ylabel('E_{SOLUTE-SOLVENT} (cm^{-1})');
xlim([-inf inf]);ylim([-inf inf]);
set(gca,'fontsize',labfont,'linewidth',lwidax);

figure;plot(r-norm(r0),EC,'b-','linewidth',4);
hold on
plot([-0.35 -0.35],[min(EC) max(EC)],'k:','linewidth',4);
% plot([-2.19 -2.19],[min(EC) max(EC)],'r:','linewidth',4);
% plot([-2.1 -2.1],[min(EC) max(EC)],'k:','linewidth',4);
plot([-0.19 -0.19],[min(EC) max(EC)],'r:','linewidth',4);
axis square;xlabel('dr (Angs)');ylabel('E_{SOLUTE-SOLVENT} (cm^{-1})');
xlim([-inf inf]);ylim([-inf inf]);
set(gca,'fontsize',labfont,'linewidth',lwidax);

figure;plot(r-norm(r0),ESS,'b-','linewidth',4);
hold on
plot([-0.35 -0.35],[min(ESS) max(ESS)],'k:','linewidth',4);
% plot([-2.19 -2.19],[min(ESS) max(ESS)],'r:','linewidth',4);
% plot([-2.1 -2.1],[min(ESS) max(ESS)],'k:','linewidth',4);
plot([-0.19 -0.19],[min(ESS) max(ESS)],'r:','linewidth',4);
axis square;xlabel('dr (Angs)');ylabel('E_{SOLUTE-SOLVENT} (cm^{-1})');
xlim([-inf inf]);ylim([-inf inf]);
set(gca,'fontsize',labfont,'linewidth',lwidax);

% scatter3(r*n(1),r*n(2),r*n(3))
% hold on
% scatter3(px(926),py(926),pz(926))
% hold on
% scatter3(hx(1:5),hy(1:5),hz(1:5))
% hold on
% scatter3(hx(48)+r*n(1),hy(48)+r*n(2),hz(48)+r*n(3))
% hold on
% scatter3(hx(49)+r*n(1),hy(49)+r*n(2),hz(49)+r*n(3))
