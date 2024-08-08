t3=0:dt/1000:150-dt/1000;
s90=SS90(12493:49992);
s100=SS100(12493:49992);
s110=SS110(12493:49992);
s120=SS120(12493:49992);
lt=length(s120);

xs90=xcorr(s90-mean(s90));xs90=xs90(lt:end)/lt;
xs100=xcorr(s100-mean(s100));xs100=xs100(lt:end)/lt;
xs110=xcorr(s110-mean(s110));xs110=xs110(lt:end)/lt;
xs120=xcorr(s120-mean(s120));xs120=xs120(lt:end)/lt;
sapod=apodfun(t3,0,0,60,60,false);

nf=2^(nextpow2(lt)+2);
fxs90=real(fftshift(fft(xs90.*sapod,nf)));
fxs100=real(fftshift(fft(xs100.*sapod,nf)));
fxs110=real(fftshift(fft(xs110.*sapod,nf)));
fxs120=real(fftshift(fft(xs120.*sapod,nf)));
w=[-1/2:1/nf:1/2-1/nf]*2*pi/dt/0.00003;dw=abs(w(1)-w(2));w=w-dw/2;

figure;plot(t3,s90,t3,s100,t3,s110,t3,s120);
figure;plot(t3,xs90/xs90(1).*sapod,t3,xs100/xs100(1).*sapod,t3,...
    xs110/xs110(1).*sapod,t3,xs120/xs120(1).*sapod,'linewidth',2);
figure;plot(w,fxs90/max(fxs90),w,fxs100/max(fxs100),w,...
    fxs110/max(fxs110),w,fxs120/max(fxs120),'linewidth',2);
xlim([0 50]);
