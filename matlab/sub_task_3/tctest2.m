clear c1 ix x y it
ix=1:500;
dt=4;
x=ix*dt-dt;x=x/1000;
x2=1:500;x2=x2/1000;
lx2=length(x2);
dt2=4;
lx=length(x);
muy=mean(mbcSS(ix));
y=mbcSS(ix);y=y-muy;
apod=apodfun(x,0,0,60,60,false);
N=muy^3;
tic
for i1=1:500
    t1=ix-i1;
    y1=circshift(y,[0,i1]);y1(t1<0)=0;
    for i2=1:500
            t2=ix-i1-i2;
%             t2=ix-i2;
            y2=circshift(y,[0,i1+i2]);y2(t2<0)=0;
%             y2=circshift(y,[0,i2]);y2(t2<0)=0;
            c1(i1,i2)=sum(y.*y1.*y2);
%             c1(i2,i1)=c1(i1,i2);
    end
end
toc
figure;contourf(x2,x2,c1/N,50,'edgecolor','none');axis square;colormap jet;colorbar;%xlim([0 100]);ylim([0 100])
fc1=fftshift(fft2(c1));fc1=fc1/max(max(fc1));
omega=((-1/2:1/lx2:1/2-1/lx2)+1/lx2/2)*2*pi/dt2/0.00003;omega=omega-abs(omega(1)-omega(2))/2;
figure;contourf(omega,omega,real(fc1),50,'edgecolor','none');axis square;
colormap jet;colorbar;%caxis([-1 1]);xlim([-10 10]);ylim([-10 10]);

y=mbcSS(1:250);y=[fliplr(y) y];y=y-mean(y);
xc=xcorr(y);xc=xc(500:999);
xc2=xcorr(y,y.^2);xc2=xc2(500:999);
xc1=xcorr(y.^2,y);xc1=xc1(500:999);
xcd=xcorr(y,y.^2);xcd=xcd(500:999);
% figure;plot(x2,c1(1,:),x2,xc1,x2,c1(:,1),x2,xc2);
figure;plot(x2,c1(1,:),x2,c1(:,1),x2,xc1,x2,diag(c1),x2,xcd);
figure;plot(x2,xc1,x2,xc2);

% figure;contourf(yy,50,'edgecolor','none');axis square;colormap jet;colorbar;
% clear yy cc tt
% xc=xcorr(y,y.^2);xc=xc(1000:1999);
%     for ii=1:1000
%             tt=ix-ii;
%             yy=circshift(y,[0,ii]);yy(tt<0)=0;
%             cc(ii)=sum(y.*yy.^2);
%     end
% figure;plot(x,cc,x,xc);
% 
% 
% xc1=xcorr(y,y.^2);xc1=xc1(1000:1999);
% xc2=xcorr(y.^2,y);xc2=xc2(1000:1999);
% figure;plot(x,xc1,x,-xc2);
