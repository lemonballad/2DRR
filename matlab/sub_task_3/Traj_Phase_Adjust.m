traj=mbcSS-mean(mbcSS);dt=(t(2)-t(1))*1000;
temp=xcorr(traj);n=length(traj);t_corr=temp(n:end)/n/dt;
omega=((-1/2:1/n:1/2-1/n)+1/n/2)*2*pi/dt/0.00003;%omega=omega-abs(omega(1)-omega(2));
apod=apodfun(t,0,0,100,100,false);
f_corr=fftshift(fft(t_corr));
efun=exp(-2*pi*5308/200*0.00003*omega);phi=atan((1-efun)./(1+efun));
figure;plot(omega,phi);
test_f=abs(f_corr).*exp(1i*phi);
fun = @(x,test_f)angle(x(1)*test_f+x(2));
x0 = [1+1i,1+1i];x = lsqcurvefit(fun,x0,test_f,phi);
test_f2=x(1)*test_f+x(2);test_f2=test_f2-mean(test_f2(m1:end));
figure;plot(omega,test_f,omega,test_f2,omega,f_corr);xlim([0 inf])
sum(omega(omega>=0).*test_f(omega>=0))/sum(real(test_f(omega>=0)))
sum(omega(omega>=0).*test_f2(omega>=0))/sum(real(test_f2(omega>=0)))
sum(omega(omega>=0).*f_corr(omega>=0))/sum(real(f_corr(omega>=0)))
figure;plot(omega,f_corr/sum(real(f_corr)),omega,test_f2/sum(real(test_f2)));

