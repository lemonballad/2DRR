clear all
labfont=16;lwidax=2;c=0.00003;
nt=128;nf=1024;nf2=1024;dt=10;
rng default;  % For reproducibility
r=2*rand(1,nt)-1;%r=r-mean(r);
x=0:dt:dt*(nt-1);
w1=500;w2=150;tau=1000;
f1=sin(c*w1*x);%
f2=exp(-x/tau).*sin(c*w2*x);%+sin(c*w2*x)+sin(c*w2*x);
f3=sin(c*w2*x);
figure;plot(x,f1,'k-',x,f2,'r-','linewidth',2);
axis square;set(gca,'fontsize',labfont,'linewidth',lwidax);
xlabel('t (fs)');ylabel('V(t) (cm^{-1})');%title(['kappa = ' num2str(kappa)]);

mu1=mean(f1);
mu2=mean(f2);
xf1=xcorr(f1-mu1);xf1=xf1(nt:end)/nt;
xf2=xcorr(f2-mu2);xf2=xf2(nt:end)/nt;
figure;plot(x,xf1,'k-',x,xf2,'r-','linewidth',2);
axis square;set(gca,'fontsize',labfont,'linewidth',lwidax);
xlabel('t (fs)');ylabel('V(t) (cm^{-1})');%title(['kappa = ' num2str(kappa)]);

for i1=1:nt
    t1=x/dt-i1;
    y1=circshift(f1,[0,i1]);y1(t1<0)=0;
    for i2=1:nt
            t2=x/dt-i1-i2;
            y2=circshift(f2,[0,i1+i2]);y2(t2<0)=0;
            for i3=1:nt
                t3=x/dt-i1-i2-i3;%x/dt-i1-i2-i3;
                y3=circshift(f3,[0,i1+i2+i3]);y3(t3<0)=0;
                c2(i1,i2,i3)=sum(f1.*y1.*y2.*y3)/nt/nt/nt;
            end
    end
end
ii=1;
c1=permute(c2(:,ii,:),[1 3 2]);
apod=apodfun(x,0,0,tau,tau,false);
apod2=repmat(apod,nt,1);apod2=apod2.*apod2';

% figure;contourf(x/1000,x/1000,c1.*apod2,50,'edgecolor','none');axis square;colormap jet;colorbar;
% set(gca,'fontsize',labfont,'linewidth',lwidax);legend('C(\tau_1,\tau_2)');
% xlabel('t_1 (ps)');ylabel('t_1+t_2 (ps)');%title(['kappa = ' num2str(kappa)]);

fc1=fftshift(fft2(c1.*apod2,nf2,nf2));
omega2=((-1/2:1/nf2:1/2-1/nf2)+1/nf2/2)*2*pi/dt/c;df2=abs(omega2(1)-omega2(2));
omega2=omega2-df2/2;
orange2=1:nf;%ceil(nf2/2+1-500/df2):ceil(nf2/2+1+500/df2);
figure;contourf(omega2(orange2),omega2(orange2),abs(fc1(orange2,orange2)),50,'edgecolor','none');
axis square;colormap jet;colorbar;xlim([-500 500]);ylim([-500 500]);
set(gca,'fontsize',labfont,'linewidth',lwidax);legend('2D SD');
xlabel('\omega_1 (cm^{-1})');ylabel('\omega_2 (cm^{-1})');
% title(['kappa = ' num2str(kappa)]);

apod=apodfun(x,0,0,tau,tau,false);
[a1,a2,a3]=ndgrid(apod,apod,apod);apod3=a1.*a2.*a3;

fct=fftshift(fftshift(fft(fft(c2.*apod3,nf2,1),nf2,3),1),3);
[~,iw1]=min(abs(omega2-w1));
fct=permute(fct(iw1,:,:),[3,2,1]);
omega2=((-1/2:1/nf2:1/2-1/nf2)+1/nf2/2)*2*pi/dt/c;df2=abs(omega2(1)-omega2(2));
omega2=omega2-df2/2;
orange2=1:nf;%ceil(nf2/2+1-500/df2):ceil(nf2/2+1+500/df2);
figure;contourf(x,omega2(orange2),abs(fct(orange2,:)),50,'edgecolor','none');
axis square;colormap jet;colorbar;xlim([0 2500]);ylim([0 500]);
set(gca,'fontsize',labfont,'linewidth',lwidax);legend('2D SD');
xlabel('\omega_1 (cm^{-1})');ylabel('\omega_2 (cm^{-1})');
% title(['kappa = ' num2str(kappa)]);
