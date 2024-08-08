nt=500;
nf=1024;
nf2=4096;
dt=4;
mux=10;
muy=10;
sx=10;
sy=10;
p=-0.8;
mu=[mux,muy];
SIGMA = [sx^2 p*sx*sy;p*sx*sy sy^2];
rng default;  % For reproducibility
r = mvnrnd(mu,SIGMA,nt);
% r=randn(2,500);r=cholcov(SIGMA)*r;
x=0:dt:dt*(nt-1);
f=r(:,1)'+r(:,2)';
tau=50;
% arf=exp(-(dt)/2/tau).*circshift(f,[0,1]);arf(1)=f(1);
% f=f+arf;
muf=mean(f);
vuf=var(f);
f=f-muf;
figure;plot(x,f);
apod=apodfun(x,0,0,500,500,false);
apod2=repmat(apod,nt,1);apod2=apod2.*apod2';
xf=xcorr(f);xf=xf(nt:(2*nt-1))/nt;
figure;plot(x,xf.*apod);
for i1=1:nt
    t1=x/dt-i1;
    y1=circshift(f,[0,i1]);y1(t1<0)=0;
    for i2=1:nt
            t2=x/dt-i1-i2;
            y2=circshift(f,[0,i1+i2]);y2(t2<0)=0;
            c1(i1,i2)=sum(f.*y1.*y2);
    end
end
c1=c1/muf^3;
figure;contourf(1:nt,1:nt,c1,50,'edgecolor','none');axis square;colormap jet;colorbar;
fc=fftshift(fft(xf.*apod,nf));fc=fc/max(fc);
fc1=fftshift(fft2(c1,nf2,nf2));fc1=fc1/max(max(fc1));
omega=((-1/2:1/nf:1/2-1/nf)+1/nf/2)*2*pi/dt/0.00003;df=abs(omega(1)-omega(2));
omega=omega-df/2;
omega2=((-1/2:1/nf2:1/2-1/nf2)+1/nf2/2)*2*pi/dt/0.00003;df2=abs(omega2(1)-omega2(2));
omega2=omega2-df2/2;
orange=ceil(nf/2+1-0/df):ceil(nf/2+1+1000/df);
orange2=ceil(nf2/2+1-500/df2):ceil(nf2/2+1+500/df2);
figure;plot(omega,abs(fc));axis square
figure;contourf(omega2(orange2),omega2(orange2),abs(fc1(orange2,orange2)),50,'edgecolor','none');axis square;
colormap jet;colorbar;%xlim([-1000 1000]);ylim([-1000 1000]);


