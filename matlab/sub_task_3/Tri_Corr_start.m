clear c1 ix x y
ix=1:1:49992;
dt=4;
x=ix*dt-dt;x=x/1000;
lx=length(x);
muy=mean(SS90(ix));
y=SS90(ix);y=y-muy;%zeros(1,lx);y([1:50 100:110])=1;%0.7*exp(-x.^2/2/1000^2)+0.3*exp(-(x).^2/2/2000^2);%exp(-x/300);
apod=apodfun(x,0,0,60,60,false);ia0=find(apod<10^-3,1);if isempty(ia0),ia0=lx;end;
N=muy^3;
ix=ix(1:ia0);
x=x(1:ia0);
dt=1000*(x(2)-x(1));
tic
for i1=1:ia0
    for i2=i1:ia0
        c1(i1,i2)=0;
        for ii=1:lx
            t1=ii-i1;
            t2=ii-i2;
            if t1<0||t2<0
                c1(i1,i2)=c1(i1,i2)+0;%
            else
                c1(i1,i2)=c1(i1,i2)+y(ii)*y(t1+1)*y(t2+1);
            end
        end
        c1(i1,i2)=c1(i1,i2)*apod(i1)*apod(i2);
        c1(i2,i1)=c1(i1,i2);
    end
end
% c1(ia0+1:lx,ia0+1:lx)=0;
toc
figure;contourf(x,x,c1/muy^3,50,'edgecolor','none');axis square;colormap jet;colorbar;xlim([0 100]);ylim([0 100])
fc1=fftshift(fft2(c1));fc1=fc1/max(max(fc1));
lx=length(x);
omega=((-1/2:1/lx:1/2-1/lx)+1/lx/2)*2*pi/dt/0.00003;%omega=omega-abs(omega(1)-omega(2));
figure;contourf(omega,omega,real(fc1),50,'edgecolor','none');axis square;
% map=[238/255 130/255 238/255;0 0 1;1 1 1;1 1 0;1 0 0];map=interp1(1:5,map,1:0.01:5);
map=[0 0 1;0 1 1;1 1 1;1 1 0;1 0 0];map=interp1(1:5,map,1:0.01:5);
colormap(map);colorbar;caxis([-1 1]);xlim([-100 100]);ylim([-100 100]);

fc2=real(fftshift(fft(c1/N/2/dt,[],1),1));fc2=fc2/max(max(fc2));
figure;contourf(x,omega,fc2,50,'edgecolor','none');axis square;colormap jet;colorbar;xlim([0 x(end)]);ylim([0 20])

if false
xc=xcorr(y);xc=xc(1001:end);

clear c1
x=0:1000;
y=mbcSS(x+1);y=y-mean(y);%exp(-x/300);
for i1=1:1001
        c1(i1)=0;
        for ii=1:1001
            t1=x(ii)-x(i1);
            if t1<0
                c1(i1)=c1(i1)+0;
                y1(ii,i1)=0;
            else
                c1(i1)=c1(i1)+y(ii)*y(ii-i1+1);
                y1(ii,i1)=y(ii-i1+1);
            end
            tt(ii,i1)=t1;
        end
    end
figure;plot(x,xc/max(xc),x,c1/max(c1));axis square;

clear x y c1
x=-100:100;
y=zeros(1,101);y(1:20)=1;y=[fliplr(y(2:end)) y];%exp(-x.^2/2/100^2);%mbcSS(x+1);y=y-mean(y);%exp(-x/300);
for i1=1:201
    t1=x(i1);
    for i2=1:201
        t2=x(i2);
        c1(i1,i2)=0;
        for ii=1:201
            t=x(ii);
            tau1=t-t1;itau1=tau1+101;
            tau2=t-t2;itau2=tau2+101;
            if abs(tau1)>100||abs(tau2)>100
                c1(i1,i2)=c1(i1,i2)+0;
            else
                c1(i1,i2)=c1(i1,i2)+y(ii)*y(itau1)*y(itau2);
            end
        end
    end
end
figure;contourf(c1,50,'edgecolor','none');axis square;colormap jet;colorbar;

end