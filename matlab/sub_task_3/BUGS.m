    labfont=16;lwidax=2;n=1;
    dt=4*n;t=0:dt/1000:(length(mbcSS)-1)*dt/1000;
    A=moviein(10000);
    hf=figure;
for ii=1:500:100001
    trange=1:n:10000;trange=trange+ii-1;
    Solu_Solv_trj=mbcSS(trange);
    Nf=length(Solu_Solv_trj);
    mean_int=mean(Solu_Solv_trj);var_int=var(Solu_Solv_trj);std_int=sqrt(var_int);
    temp=xcorr(Solu_Solv_trj-mean_int,'coeff');t_corr=temp(Nf:end);
    time=t(trange)-t(trange(1));f=fit(time',t_corr','exp2','StartPoint',[0.2 -4 0.8 -12]);
    apod=apodfun(time,0,0,1,1,false);%figure;plot(time,t_corr,time,t_corr.*apod,time,apod);
    t_corr=t_corr.*apod;t_corr(apod==0)=[];t_corr=[t_corr 0];apod(apod==0)=[];apod=[apod 0];
    f_corr=abs(fftshift(fft(t_corr.*apod)));f_corr=f_corr-f_corr(end);
    Nf=length(t_corr);omega=((-1/2:1/Nf:1/2-1/Nf)+1/Nf/2)*2*pi/dt/0.00003;
    plot(omega,f_corr/max(f_corr),'linewidth',lwidax);
    xlim([0 2000]);ylim([0 inf]);xlabel('\omega (cm^{-1})');ylabel('C(\omega)');title(['Carboxy-Myoglobin: t=' num2str(trange(1)) '-' num2str(trange(end))]);
    set(gca,'fontsize',labfont,'linewidth',lwidax);axis square;
    A(:,ii)=getframe(hf);
end    
