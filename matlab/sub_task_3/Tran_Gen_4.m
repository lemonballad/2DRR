clear all
labfont=16;
lwidax=2;
nt=1024;
nf=2048;
nf2=4096;
dt=10;
c=0.00003;
rng default;  % For reproducibility
r=rand(2,nt);%r=r-mean(r);
x=0:dt:dt*(nt-1);
tau=10^3;
kappa=1;
fluc=((1/tau/c)/kappa);
f=r(:,:);
p=-0.8;SIGMA=[1 p;p 1];
f=(cholcov(SIGMA)*r);
arf=zeros(2,nt);
for ii=0:nt
    arf=arf+...
        exp(-kappa^-2*exp(-(ii*dt)/tau)).*...
        circshift(f,[0,ii]);
end
f=(arf);
muf=mean(f');
f=sum(f)-sum(muf);
figure;plot(x,f,'k-','linewidth',2);
axis square;set(gca,'fontsize',labfont,'linewidth',lwidax);
xlabel('t (fs)');ylabel('V(t) (cm^{-1})');title(['kappa = ' num2str(kappa)]);

apod=apodfun(x,0,0,1/fluc/c,1/fluc/c,false);
apod2=repmat(apod,nt,1);apod2=apod2.*apod2';
xf=xcorr(f);xf=xf(nt:(2*nt-1))/nt.*apod;
figure;plot(x,xf,'k-',[tau tau],[min(xf) max(xf)],'r:','linewidth',2);
axis square;set(gca,'fontsize',labfont,'linewidth',lwidax);
xlabel('t (fs)');ylabel('C(t) (cm^{-2})');title(['\kappa = ' num2str(kappa)]);
legend(['\Delta^2 = ' num2str(round(fluc^2)) ' cm^{-2}'],['\tau = ' num2str(tau) ' fs']);
ylim([min(xf) max(xf)]);xlim([0 10*tau]);

w_vib=0;
omega=((-1/2:1/nf:1/2-1/nf)+1/nf/2)*2*pi/dt/c;df=abs(omega(1)-omega(2));
omega=omega-df/2;
tic
for i1=1:nt
    t1=x/dt-i1+1;
    y1=circshift(f,[0,i1-1]);y1(t1<0)=0;
    for i2=1:nt
            t2=x/dt-i1-i2+2;
            y2=circshift(f,[0,i1+i2-2]);y2(t2<0)=0;
            c1(i1,i2)=sum(f.*y1.*y2);
    end
end
toc
c1=c1/nt/nt;

figure;contourf(x/1000,x/1000,c1,50,'edgecolor','none');axis square;colormap jet;colorbar;
set(gca,'fontsize',labfont,'linewidth',lwidax);legend('C(\tau_1,\tau_2)');
xlabel('t_1 (ps)');ylabel('t_1+t_2 (ps)');title(['kappa = ' num2str(kappa)]);

fc=fftshift(fft(xf/5308/5308,nf));
fc1=fftshift(fft2(c1.*apod2,nf2,nf2));
orange=ceil(nf/2+1-0/df):ceil(nf/2+1+1000/df);
omega2=((-1/2:1/nf2:1/2-1/nf2)+1/nf2/2)*2*pi/dt/c;df2=abs(omega2(1)-omega2(2));
omega2=omega2-df2/2;
orange2=ceil(nf2/2+1-500/df2):ceil(nf2/2+1+500/df2);
figure;contourf(omega2(orange2),omega2(orange2),abs(fc1(orange2,orange2)),50,'edgecolor','none');
axis square;colormap jet;colorbar;xlim([-500 500]);ylim([-500 500]);
set(gca,'fontsize',labfont,'linewidth',lwidax);legend('2D SD');
xlabel('\omega_1 (cm^{-1})');ylabel('\omega_2 (cm^{-1})');
title(['kappa = ' num2str(kappa)]);

w=((-1/2:1/nf:1/2-1/nf)+1/nf/2)/dt/c;df=abs(w(1)-w(2));
w=w-df/2;
fc=fc-fc(end);
figure;plot(omega,real(fc),'k-'...
    ,[fluc fluc],[0 max(real(fc))],'r:','linewidth',2);
axis square;xlim([0 500]);set(gca,'fontsize',labfont,'linewidth',lwidax);
xlabel('\omega (cm^{-1})');legend('SD',['\Delta = ' num2str(round(fluc)) ' cm^{-1}']);
title(['kappa = ' num2str(kappa)]);ylim([0 max(real(fc))]);

