clear c1 ix x y it
ix=1:5:50000;
dt=4;
x=ix*dt-dt;x=x/1000;
lx=length(x);
muy=mean(mbdSS(ix));
y=mbdSS(ix);y=y-muy;
apod=apodfun(x,0,0,60,60,false);ia0=find(apod<10^-3,1);if isempty(ia0),ia0=lx;end;
N=muy^3;
c1=zeros(ia0,ia0,'double');
% cy=zeros(ia0,lx,'double');
% for ii=1:ia0
%     cy(ii,:)=circshift(y,0,ii);cy(ii,ix-ii<0)=0;
% end
tic
for i1=1:ia0
    t1=ix-i1;
    y1=circshift(y,[0,i1]);y1(t1<0)=0;%y1=y1(it);%cy(i1,:);
    for i2=i1:ia0
            t2=ix-i2;%it=t1>=0&t2>=0;
            y2=circshift(y,[0,i2]);y2(t2<0)=0;%y2=y2(it);%cy(i2,:);
            c1(i1,i2)=sum(y.*y1.*y2);
%             c1(i1,i2)=sum(y.*cy(i1,:).*cy(i2,:));
        c1(i1,i2)=c1(i1,i2)*apod(i1)*apod(i2);
        c1(i2,i1)=c1(i1,i2);
    end
end
% c1(ia0+1:lx,ia0+1:lx)=0;
toc
% ix=ix(1:ia0);
% x=x(1:ia0);
% dt=1000*(x(2)-x(1));
figure;contourf(x,x,c1/N,50,'edgecolor','none');axis square;colormap jet;colorbar;xlim([0 100]);ylim([0 100])
fc1=fftshift(fft2(c1));fc1=fc1/max(max(fc1));
lx=length(x);
omega=((-1/2:1/lx:1/2-1/lx)+1/lx/2)*2*pi/dt/0.00003;omega=omega-abs(omega(1)-omega(2))/2;
fc1=fc1(omega>-100&omega<100,omega>-100&omega<100);
fc2=real(fftshift(fft(c1/N/2/dt,[],1),1));fc2=fc2/max(max(fc2));fc2=fc2(omega>-100&omega<100,:);
omega=omega(omega>-100&omega<100);
figure;contourf(omega,omega,real(fc1),50,'edgecolor','none');axis square;
% map=[238/255 130/255 238/255;0 0 1;1 1 1;1 1 0;1 0 0];map=interp1(1:5,map,1:0.01:5);
map=[0 0 1;0 1 1;1 1 1;1 1 0;1 0 0];map=interp1(1:5,map,1:0.01:5);
colormap(map);colorbar;caxis([-1 1]);xlim([-100 100]);ylim([-100 100]);

figure;contourf(x,omega,fc2,50,'edgecolor','none');axis square;colormap jet;colorbar;xlim([0 x(end)]);ylim([0 200])

if false
xc=xcorr(y);xc=xc(1001:end);

clear c1
x=0:1000;
y=mbcSS(x+1);y=y-mean(y);%exp(-x/300);
for i1=1:1001
        c1(i1)=0;
            t1=x-x(i1);y1=circshift(y,[0,i1]);y1(t1<0)=0;
            yy(:,i1)=y1;
                c1(i1)=sum(y.*y1);
            ttt(:,i1)=t1;
end
figure;plot(x,xc/max(xc),x,c1/max(c1));axis square;



yz=y((t1>=0)&(t2>=0));
end