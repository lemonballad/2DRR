clear all
labfont=16;
lwidax=2;
nt=2000;
nf=2048;
nf2=4096;
dt=10;
c=0.00003;
x=0:dt:dt*(nt-1);
tau=10^3;
kappa=1.5;
fluc=((1/tau/c)/kappa);
apod=apodfun(x,0,0,700,700,false);
f=exp(-kappa^-2*(exp(-x/tau)+x/tau-1));
figure;plot(x,f,'k-',[tau tau],[min(f) max(f)],'r:','linewidth',2);
axis square;set(gca,'fontsize',labfont,'linewidth',lwidax);
xlabel('t (fs)');ylabel('C(t) (cm^{-2})');title(['\kappa = ' num2str(kappa)]);
legend(['\Delta^2 = ' num2str(round(fluc^2)) ' cm^{-2}'],['\tau = ' num2str(tau) ' fs']);
ylim([min(f) max(f)]);xlim([0 10*tau]);

xf=xcorr(f);xf=xf(nt:end)/nt;
figure;plot(x,xf,'k-',[tau tau],[min(xf) max(xf)],'r:','linewidth',2);
axis square;set(gca,'fontsize',labfont,'linewidth',lwidax);
xlabel('t (fs)');ylabel('C(t) (cm^{-2})');title(['\kappa = ' num2str(kappa)]);
legend(['\Delta^2 = ' num2str(round(fluc^2)) ' cm^{-2}'],['\tau = ' num2str(tau) ' fs']);
ylim([min(xf) max(xf)]);
% g=cumtrapz(x,cumtrapz(x,f/5308/5308));
% figure;plot(g);
% w_vib=0;
% w=((-1/2:1/nf:1/2-1/nf)+1/nf/2)/dt/c;df=abs(w(1)-w(2));
% w=w-df/2;
omega=((-1/2:1/nf:1/2-1/nf)+1/nf/2)*2*pi/dt/c;df=abs(omega(1)-omega(2));
omega=omega-df/2;
% sig_t=exp(-1i*2*pi*w_vib*c*x-g);
% sig_w=fftshift(fft(sig_t,nf));sig_w=sig_w-sig_w(end);
% figure;plot(w,real(sig_w)/max(real(sig_w)),'linewidth',2);
% set(gca,'fontsize',labfont,'linewidth',lwidax);

fc=fftshift(fft(f,nf));fc=fc-fc(end);%fc=fc/max(real(fc));
% rsw=max(real(fc))*real(sig_w)/max(real(sig_w));
% figure;plot(w,rsw,'b-',omega,real(fc),'k-'...
figure;plot(omega,real(fc),'k-'...
    ,[fluc fluc],[0 max(real(fc))],'r:','linewidth',2);
axis square;xlim([0 100]);set(gca,'fontsize',labfont,'linewidth',lwidax);
xlabel('\omega (cm^{-1})');legend('SD',['\Delta = ' num2str(round(fluc)) ' cm^{-1}']);
% legend('\sigma','SD',['\Delta = ' num2str(round(fluc)) ' cm^{-1}']);
% title(['kappa = ' num2str(kappa)]);ylim([0 max([max(rsw) max(real(fc))])]);
title(['kappa = ' num2str(kappa)]);ylim([0 max(real(fc))]);

ff=f;%fluc*sqrt(2)*exp(-x/tau);
for i1=1:nt
    t1=x/dt-i1;
    y1=circshift(ff,[0,i1]);y1(t1<0)=0;
    for i2=1:nt
            t2=x/dt-i1-i2;
            y2=circshift(ff,[0,i1+i2]);y2(t2<0)=0;
            c1(i1,i2)=sum(ff.*y1.*y2);
    end
end
c1=c1/nt/nt;

omega2=((-1/2:1/nf2:1/2-1/nf2)+1/nf2/2)*2*pi/dt/c;df2=abs(omega2(1)-omega2(2));
omega2=omega2-df2/2;
orange2=ceil(nf2/2+1-1000/df2):ceil(nf2/2+1+1000/df2);

figure;contourf(x/1000,x/1000,c1,50,'edgecolor','none');axis square;colormap jet;colorbar;
set(gca,'fontsize',labfont,'linewidth',lwidax);legend('C(\tau_1,\tau_2)');
xlabel('t_1 (ps)');ylabel('t_1+t_2 (ps)');title(['kappa = ' num2str(kappa)]);
xlim([0 5]);ylim([0 5]);

fc1=fftshift(fft2(c1,nf2,nf2));%fc1=fc1/max(max(real(fc1)));
orange=ceil(nf/2+1-0/df):ceil(nf/2+1+1000/df);
figure;contourf(omega2(orange2),omega2(orange2),abs(fc1(orange2,orange2)),50,'edgecolor','none');
axis square;colormap jet;colorbar;xlim([-500 500]);ylim([-500 500]);
set(gca,'fontsize',labfont,'linewidth',lwidax);legend('2D SD');
xlabel('\omega_1 (cm^{-1})');ylabel('\omega_2 (cm^{-1})');
title(['kappa = ' num2str(kappa)]);


