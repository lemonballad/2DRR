xc=zeros(20,2500,'double');
mu_avg=zeros(20,'double');
f_avg=zeros(20,'double');
for ii=1:47475
    range=(1:2500)+(ii-1);
    for jj=2:20
        temp=wvn(jj,range);
        mu=mean(temp);
        mu_avg(jj)=mu_avg(jj)+mu;
        muf=mean(f(jj,range),2);
        f_avg(jj)=f_avg(jj)+muf;
        temp=xcorr(temp-mu);
        temp2=temp(2500:end)/2500;
        xc(jj,:)=xc(jj,:)+temp2;
    end
    if mod(ii,1000)==0
        ii
    end;
end
f_avg=f_avg/47475;
mu_avg=mu_avg/47475;
xc=xc/47475;

nf=4096;dt=20;c=0.00003;
omega=((-1/2:1/nf:1/2-1/nf)+1/nf/2)*2*pi/dt/c;df=abs(omega(1)-omega(2));omega=omega-df/2;

t=0:0.02:2500*0.02-0.02;
apod=repmat(apodfun(t,0,0,5,5,false),20,1);

fxc=real(fftshift(fft(xc.*apod,nf,2)));
w=30001:50000;
spec=zeros(1,20000);
for ii=[2:10 12:20]
    w_temp=omega+mu_avg(ii);
    sfxc=fxc(ii,:);%sfxc=sfxc/max(sfxc);%*f_avg(ii);
    specspec+interp1(w_temp,smooth(sfxc),w,'spline',0);
end
