apod=apodfun(t,0,0,60,60,false);n=length(mbcSS);
dt=4;omega=((-1/2:1/n:1/2-1/n)+1/n/2)*2*pi/dt/0.00003;omega=omega-abs(omega(1)-omega(2));
traj=mbcSS(indices)-mean(mbcSS);temp=xcorr(traj);ntraj=length(traj);t_corr_C=temp(ntraj:end)/ntraj/4;
traj=mbdSS(indices)-mean(mbdSS);temp=xcorr(traj);ntraj=length(traj);t_corr_D=temp(ntraj:end)/ntraj/4;
% [~,m1]=min(abs(omega-2000));[~,m2]=min(abs(omega-3000));
% [~,m]=min(abs(omega-200));[~,m0]=min(abs(omega-0));lo=length(omega(m1:m2));
f_corr_C=real(fftshift(fft(t_corr_C.*apod)));
f_corr_C=f_corr_C-f_corr_C(end);%sum(f_corr_C(m1:m2))/lo;
f_corr_C=f_corr_C/max(f_corr_C);
f_corr_D=real(fftshift(fft(t_corr_D.*apod)));
f_corr_D=f_corr_D-f_corr_D(end);%-sum(f_corr_D(m1:m2))/lo;
f_corr_D=f_corr_D/max(f_corr_D);

CMass_C=dlmread('CMass_mbc.txt');paramEsts_C=dlmread('params_mbc.txt');
CMass_D=dlmread('CMass_mbd.txt');paramEsts_D=dlmread('params_mbd.txt');

pdf_C=@(x,p,mu1,mu2,sigma1,sigma2)p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);                 
pdf_C=pdf_C(omega,paramEsts_C(1),paramEsts_C(2),paramEsts_C(3),paramEsts_C(4),paramEsts_C(5));
dist_D=fitdist(CMass_D,'Lognormal');pdf_D=pdf(dist_D,omega);

dlmwrite('time.txt',t);dlmwrite('omega.txt',omega);
dlmwrite('mbdSS.txt',mbdSS);dlmwrite('f_corr_D.txt',f_corr_D);
dlmwrite('CMass_D.txt',CMass_D);dlmwrite('pdf_D.txt',pdf_D);
dlmwrite('mbcS.txt',mbcSS);dlmwrite('f_corr_C.txt',f_corr_C);
dlmwrite('CMass_C.txt',CMass_C);dlmwrite('pdf_C.txt',pdf_C);