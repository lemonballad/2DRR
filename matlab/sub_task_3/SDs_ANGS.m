SS90=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an90\E_SS\SS_NOCO_90.txt');
SS95=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an95\E_SS\SS_NOCO_95.txt');
SS100=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an100\E_SS\SS_NOCO_100.txt');
SS105=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an105\E_SS\SS_NOCO_105.txt');
SS110=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an110\E_SS\SS_NOCO_110.txt');
SS120=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an120\E_SS\SS_NOCO_120.txt');

SS90=SS90';SS95=SS95';SS100=SS100';SS105=SS105';SS110=SS110';SS120=SS120';
SS90=SS90(:);SS95=SS95(:);SS100=SS100(:);SS105=SS105(:);SS110=SS110(:);SS120=SS120(:);
SS90=SS90';SS95=SS95';SS100=SS100';SS105=SS105';SS110=SS110';SS120=SS120';
SS90(isnan(SS90))=[];SS95(isnan(SS95))=[];SS100(isnan(SS100))=[];SS105(isnan(SS105))=[];
SS110(isnan(SS110))=[];SS120(isnan(SS120))=[];
lt=99994;
SS90=SS90(1:lt);SS95=SS95(1+10^5:lt+10^5);SS100=SS100(1:lt);SS105=SS105(1+10^5:lt+10^5);
SS110=SS110(1:lt);SS120=SS120(1:lt);
mu90=mean(SS90);mu95=mean(SS95);mu100=mean(SS100);
mu105=mean(SS105);mu110=mean(SS110);mu120=mean(SS120);
dt=4;
t=0:dt/1000:lt*dt/1000-dt/1000;
apod=apodfun(t,0,0,60,60,false);

x90=xcorr(SS90-mu90);x90=x90(lt:end)/lt;
x95=xcorr(SS95-mu95);x95=x95(lt:end)/lt;
x100=xcorr(SS100-mu100);x100=x100(lt:end)/lt;
x105=xcorr(SS105-mu105);x105=x105(lt:end)/lt;
x110=xcorr(SS110-mu110);x110=x110(lt:end)/lt;
x120=xcorr(SS120-mu120);x120=x120(lt:end)/lt;

figure;plot(t,SS90,t,SS95,t,SS100,t,SS105,t,SS110,t,SS120,'linewidth',2);
xlim([0 400]);legend('90','95','100','105','110','120');
figure;plot(t,x90/max(x90),t,x95/max(x95),t,x100/max(x100),...
    t,x105/max(x105),t,x110/max(x110),t,x120/max(x120),'linewidth',2);
xlabel('t (ps)');ylabel('C_{SS}(t) cm^{-2}');
xlim([0 200]);legend('90','95','100','105','110','120');
figure;plot(t,x90/max(x90).*apod,t,x95/max(x95).*apod,...
t,x105/max(x105).*apod,t,x100/max(x100).*apod,...
t,x110/max(x110).*apod,t,x120/max(x120).*apod,'linewidth',2);
xlabel('t (ps)');ylabel('E_{SS}(t) cm^{-1}');
xlim([0 200]);legend('90','95','100','105','110','120');

nf=2^(nextpow2(lt)+2);
fx90=real(fftshift(fft(x90.*apod,nf)));
fx95=real(fftshift(fft(x95.*apod,nf)));
fx100=real(fftshift(fft(x100.*apod,nf)));
fx105=real(fftshift(fft(x105.*apod,nf)));
fx110=real(fftshift(fft(x110.*apod,nf)));
fx120=real(fftshift(fft(x120.*apod,nf)));
w=[-1/2:1/nf:1/2-1/nf]*2*pi/dt/0.00003;dw=abs(w(1)-w(2));w=w-dw/2;

figure;plot(w,fx90/max(fx90),w,fx95/max(fx95),w,fx100/max(fx100),...
    w,fx105/max(fx105),w,fx110/max(fx110),w,fx120/max(fx120),'linewidth',2);
xlabel('\omega (cm^{-1})');ylabel('C_{SS}(\omega) cm^{-2}');
xlim([0 50]);legend('90','95','100','105','110','120');

figure;plot(w,fx90,w,fx95,w,fx100,w,fx105,w,fx110,w,fx120,'linewidth',2);
xlabel('\omega (cm^{-1})');ylabel('C_{SS}(\omega) cm^{-2}');
xlim([0 50]);legend('90','95','100','105','110','120');



