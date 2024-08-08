nf=length(mbcSS);
nt=nf;
c=0.00003;
dt=4;
t=0:dt:nt*dt-dt;t=t/1000;
w=((-1/2:1/nf:1/2-1/nf)+1/nf/2)*2*pi/dt/c;df=abs(w(1)-w(2));
w=w-df/2;
apod=apodfun(w,0,-1/10^5/0.00003,1/10^5/c,0,true);
% fw=(1-exp(-(abs(w)-8300).^2/2/2000^2)-exp(-w.^2/2/(1)^2));fw(abs(w)>8300)=0;
fw=(1-exp(-w.^2/2/(1)^2));
plot(w,fw);
plot(t,real(ifft(ifftshift(apod))));


Cnew=ifft(ifftshift(fftshift(fft(mbcSS)).*fw));
Dnew=ifft(ifftshift(fftshift(fft(mbdSS)).*fw));

xc=xcorr(Cnew-mean(Cnew));xc=xc(nt:end)/nt;
xd=xcorr(Dnew-mean(Dnew));xd=xd(nt:end)/nt;

apod=apodfun(t,0,0,33,33,false);

fc=fftshift(fft(xc.*apod));fc=fc-fc(end);
fd=fftshift(fft(xd.*apod));fd=fd-fd(end);

lwidax=2;labfont=16;

figure;plot(t,xc*xd(1)/xc(1).*apod,'r-',t,xd.*apod,'b-','linewidth',2);xlim([0 100]);
set(gca,'fontsize',labfont,'linewidth',lwidax);legend('C(\tau_1,\tau_2)');
xlabel('t_1 (ps)');ylabel('C(t) (cm^{-2})');legend('COMb','DeoxyMb');

figure;plot(w,fc,'r-',w,fd,'b-','linewidth',2);xlim([0 50]);
set(gca,'fontsize',labfont,'linewidth',lwidax);legend('C(\tau_1,\tau_2)');
xlabel('\omega (cm^{-1})');ylabel('C(\omega)');legend('COMb','DeoxyMb');
