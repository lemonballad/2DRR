w=((-1/2:1/5001:1/2-1/5001)+1/5001/2)*2*pi/10/0.00003;df=abs(w(1)-w(2));
w=w-df/2;
apod=apodfun(t,0,0,tau,tau,false);

plot(t,x1/max(x1))

nt=200;
xc=zeros(2,nt,'double');
mu_avg=[0;0];
for ii=1:5002-nt
    range=(1:nt)+(ii-1);
    temp1=i1(range)';mu1=mean(temp1);
    temp2=i2(range)';mu2=mean(temp2);
    mu_avg=mu_avg+[mu1;mu2];
    temp1=xcorr(temp1-mu1);temp1=temp1(nt:end)/nt;
    temp2=xcorr(temp2-mu2);temp2=temp2(nt:end)/nt;
    xc(1,:)=xc(1,:)+temp1;xc(2,:)=xc(2,:)+temp2;
    if mod(ii,100)==0
        ii
    end
end
xc=xc/ii;
time=0:10:nt*10-10;
figure;plot(time,xc(1,:)/max(xc(1,:)),time,xc(2,:)/max(xc(2,:)));axis square
figure;scatter(i1,i2,'k.');axis square;xlim([3.8 4.8]);ylim([3.8 4.8]);
plot(time,sum(xc,1)/2);
w=((-1/2:1/5001:1/2-1/5001)+1/5001/2)*2*pi/10/0.00003;df=abs(w(1)-w(2));
w=w-df/2;apod=apodfun(time,0,0,250,250,false);
fc=real(fftshift(fft(xc.*[apod;apod],5001,2),2));fc(1,:)=fc(1,:)-fc(1,end);fc(2,:)=fc(2,:)-fc(2,end);
figure;plot(w,fc(1,:)/max(fc(1,:)),w,fc(2,:)/max(fc(2,:)));axis square;xlim([0 5000]);