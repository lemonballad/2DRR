for angle=[90 95 100 105 110 120];
c3=dlmread(['C:/Users/Thomas/Desktop/MBC_FE/an'...
    num2str(angle) '/E_SS/C3_' num2str(angle) '.txt']);

[nt,~]=size(c3);nf=2^14;
dt=20;
t=0:dt/1000:dt/1000*(nt-1);

% t1=repmat(t,nt,1);
% tau=t1';
% t2=tau-t1;
% alpha=interp2(t1,tau,c90,t1,t2);alpha(isnan(alpha))=0;
% alpha=alpha+alpha'-diag(diag(alpha));

labfont=16;lwidax=2;
f=((-1/2:1/nf:1/2-1/nf)+1/nf/2)*2*pi/dt/0.00003;df=abs(f(1)-f(2));f=f-df/2;
index=1:nf;
trange=1:5:3*nt/4;%ceil(100000/dt);
frange=nf/2+1:ceil(nf/2+1+20/df);
frange1=floor(nf/2+1-0/df):ceil(nf/2+1+10/df);
frange2=floor(nf/2+1-10/df):ceil(nf/2+1+10/df);
apod=apodfun(t,0,0,60,60,false);[apod,~]=ndgrid(apod,apod);apod=apod.*apod';
f1c90=fftshift(fft(c3.*apod,nf,1),1);f1c90=f1c90/max(max(f1c90));
f2c90=fftshift(fft2(c3.*apod,nf,nf));f2c90=f2c90/max(max(f2c90));

% Plot C(t1,t2)
figure;contourf(t(trange),t(trange),c3(trange,trange).*apod(trange,trange),50,'edgecolor','none');
axis square;colormap jet;colorbar;xlabel('t_1 (ps)');ylabel('t_2 (ps)');
xlim([0 t(trange(end))]);ylim([0 t(trange(end))]);
set(gca,'fontsize',labfont,'linewidth',lwidax);

% Plot C(t1,w2)
figure;contourf(t(trange),f(frange),real(f1c90(frange,trange)),100,'edgecolor','none');
axis square;colormap jet;colorbar;xlabel('t_1 (ps)');ylabel('\omega_2 (cm^{-1})');
xlim([t(min(trange)) t(max(trange))]);ylim([f(min(frange)) f(max(frange))]);
set(gca,'fontsize',labfont,'linewidth',lwidax);
% % Plot |C(t1,w2)|
% figure;contourf(t(trange),f(frange),abs(f1c90(frange,trange)),100,'edgecolor','none');
% axis square;colormap jet;colorbar;xlabel('t_1 (ps)');ylabel('\omega_2 (cm^{-1})');
% xlim([t(min(trange)) t(max(trange))]);ylim([f(min(frange)) f(max(frange))]);
% set(gca,'fontsize',labfont,'linewidth',lwidax);

% % Plot C(w1>0,w2) 
% figure;contourf(f(frange1),f(frange2),real(f2c90(frange2,frange1)),100,'edgecolor','none');
% axis square;colormap jet;colorbar;xlabel('\omega_1 (cm^{-1})');ylabel('\omega_2 (cm^{-1})');
% xlim([f(min(frange1)) f(max(frange1))]);ylim([f(min(frange2)) f(max(frange2))]);
% set(gca,'fontsize',labfont,'linewidth',lwidax);
% % Plot |C(w1>0,w2)|
% figure;contourf(f(frange1),f(frange2),abs(f2c90(frange2,frange1)),100,'edgecolor','none');
% axis square;colormap jet;colorbar;xlabel('\omega_1 (cm^{-1})');ylabel('\omega_2 (cm^{-1})');
% xlim([f(min(frange1)) f(max(frange1))]);ylim([f(min(frange2)) f(max(frange2))]);
% set(gca,'fontsize',labfont,'linewidth',lwidax);

% Plot C(w1,w2)
figure;contourf(f(frange2),f(frange2),real(f2c90(frange2,frange2)),100,'edgecolor','none');
axis square;colormap jet;colorbar;xlabel('\omega_1 (cm^{-1})');ylabel('\omega_2 (cm^{-1})');
xlim([f(min(frange2)) f(max(frange2))]);ylim([f(frange2(1)) f(frange2(end))]);
set(gca,'fontsize',labfont,'linewidth',lwidax);
% Plot |C(w1,w2)|
figure;contourf(f(frange2),f(frange2),abs(f2c90(frange2,frange2)),100,'edgecolor','none');
axis square;colormap jet;colorbar;xlabel('\omega_1 (cm^{-1})');ylabel('\omega_2 (cm^{-1})');
xlim([f(min(frange2)) f(max(frange2))]);ylim([f(frange2(1)) f(frange2(end))]);
set(gca,'fontsize',labfont,'linewidth',lwidax);

end


