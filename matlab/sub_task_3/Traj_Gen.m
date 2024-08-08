nt=1000;
nf=4096;
dt=10;
mux=10;
muy=10;
sx=5;
sy=5;
p=1;
x=randi([1 nt],1,nt)/dt;
y=randi([1 nt],1,nt)/dt;
[X,Y]=meshgrid(x,y);
eX=(x-mux)/sx;
eY=(y-muy)/sy;
% DIS=exp(-(eX.^2+eY.^2-2*p*eX.*eY));
% DIS=DIS/sum(sum(DIS));
DIS=(eX.^2+eY.^2-2*p*eX.*eY);
f=DIS;
f=f-mean(f);
% f=DIS-mean(DIS);
figure;plot(1:nt,f);
xf=xcorr(f-mean(f));xf=xf(nt:(2*nt-1))/nt;
figure;plot(1:nt,xf);

for i1=1:nt
    t1=x-i1;
    y1=circshift(f,[0,i1]);y1(t1<0)=0;
    for i2=1:nt
            t2=x-i1-i2;
            y2=circshift(f,[0,i1+i2]);y2(t2<0)=0;
            c1(i1,i2)=sum(f.*y1.*y2);
    end
end
figure;contourf(1:nt,1:nt,c1,50,'edgecolor','none');axis square;colormap jet;colorbar;
fc1=fftshift(fft2(c1,nf,nf));fc1=fc1/max(max(fc1));
omega=((-1/2:1/nf:1/2-1/nf)+1/nf/2)*2*pi/dt/0.00003;df=abs(omega(1)-omega(2));
omega=omega-df/2;
orange=ceil(nf/2+1-2000/df):ceil(nf/2+1+2000/df);
figure;contourf(omega(orange),omega(orange),real(fc1(orange,orange)),50,'edgecolor','none');axis square;
colormap jet;colorbar;%xlim([-1000 1000]);ylim([-1000 1000]);


