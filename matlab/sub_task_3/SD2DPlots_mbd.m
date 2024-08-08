betad=dlmread('C:/Users/Thomas/Desktop/SDs/mbd/mbd_C3.txt',',');
labfont=16;lwidax=2;
[nt,nf]=size(betad);dt=20;nf=2^14; %#ok<ASGLU>
f=((-1/2:1/nf:1/2-1/nf)+1/nf/2)*2*pi/dt/0.00003;df=abs(f(1)-f(2));f=f-df/2;
index=1:nf;
t=0:dt/1000:dt/1000*(nt-1);
frange=nf/2+1:ceil(nf/2+1+20/df);
trange=1:ceil(200000/dt);
frange1=floor(nf/2+1-0/df):ceil(nf/2+1+10/df);
frange2=floor(nf/2+1-10/df):ceil(nf/2+1+10/df);
apod=apodfun(t,0,0,60,60,false);[apod,~]=ndgrid(apod,apod);apod=apod.*apod';
fbetayd=fftshift(fft(betad.*apod,nf,1),1);fbetayd=fbetayd/max(max(fbetayd));
fbetad=fftshift(fft2(betad.*apod,nf,nf));fbetad=fbetad/max(max(fbetad));

figure;contourf(t(trange),t(trange),betad(trange,trange).*apod(trange,trange),50,'edgecolor','none');
axis square;colormap jet;colorbar;xlabel('t_1 (ps)');ylabel('t_1+t_2 (ps)');
xlim([0 nt*dt/1000]);ylim([0 nt*dt/1000]);
set(gca,'fontsize',labfont,'linewidth',lwidax);

figure;contourf(t(trange),f(frange),real(fbetayd(frange,trange)),100,'edgecolor','none');
axis square;colormap jet;colorbar;xlabel('t_1 (ps)');ylabel('\omega_2 (cm^{-1})');
xlim([t(min(trange)) t(max(trange))]);ylim([f(min(frange)) f(max(frange))]);
set(gca,'fontsize',labfont,'linewidth',lwidax);

figure;contourf(f(frange1),f(frange1),real(fbetad(frange1,frange1)),100,'edgecolor','none');
axis square;colormap jet;colorbar;xlabel('\omega_1 (cm^{-1})');ylabel('\omega_2 (cm^{-1})');
xlim([f(min(frange1)) f(max(frange1))]);ylim([f(min(frange1)) f(max(frange1))]);
set(gca,'fontsize',labfont,'linewidth',lwidax);

figure;contourf(f(frange2),f(frange2),real(fbetad(frange2,frange2)),100,'edgecolor','none');
axis square;colormap jet;colorbar;xlabel('\omega_1 (cm^{-1})');ylabel('\omega_2 (cm^{-1})');
xlim([f(min(frange2)) f(max(frange2))]);ylim([f(frange2(1)) f(frange2(end))]);
set(gca,'fontsize',labfont,'linewidth',lwidax);

figure;contourf(f(frange2),f(frange2),abs(fbetad(frange2,frange2)),100,'edgecolor','none');
axis square;colormap jet;colorbar;xlabel('\omega_1 (cm^{-1})');ylabel('\omega_2 (cm^{-1})');
xlim([f(min(frange2)) f(max(frange2))]);ylim([f(frange2(1)) f(frange2(end))]);
set(gca,'fontsize',labfont,'linewidth',lwidax);

y=mbdSS(1:50000);muy=mean(y);y=y-muy;
x=(1:50000)*0.004-0.004;

xc1=xcorr(y,y.^2);xc1=xc1(50000:99999)/muy^3;
xc2=xcorr(y.^2,y);xc2=xc2(50000:99999)/muy^3;
figure;plot(x(1:5:50000),(xc1(1:5:50000)),t,betad(1,:));
figure;plot(x(1:5:50000),(xc2(1:5:50000)),t,betad(:,1));