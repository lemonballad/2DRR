
n=length(mbcSS);    nstep=200;  dt=4; WindWid=100000; dw=n/nstep;
ntraj=WindWid/dt;  trange=(0:ntraj-1)*dt/1000;    srange=(0:nstep-1)*dw*dt/1000;
T_CORR=zeros(nstep,ntraj,'double'); K=zeros(1,nstep,'double');
apod=apodfun(trange,0,0,100,100,false);
for ii=1:nstep
    in=(ii-1)*dw;
    indices=1+in:ntraj+in;  indices(indices>n)=indices(indices>n)-n;
    traj=mbcSS(indices)-mean(mbcSS(indices));var_traj=var(traj);
    temp=xcorr(traj);  t_corr=temp(ntraj:end)/ntraj/var_traj;
    T_CORR(ii,1:ntraj)=t_corr.*apod;
    K(ii)=trapz(t_corr)*dt/5308/5308;
end
figure;plot(trange,T_CORR(1,:))
% figure;plot(omega,F_CORR(1,:))
TAU=1./K;
omega=((-1/2:1/ntraj:1/2-1/ntraj)+1/ntraj/2)*2*pi/dt/0.00003;omega=omega-abs(omega(1)-omega(2));
OMEGA=repmat(omega,nstep,1);
[~,m1]=min(abs(omega-200));[~,m2]=min(abs(omega-3000));
[~,m]=min(abs(omega-2000));[~,m0]=min(abs(omega-0));lo=length(omega(m1:m2));
F_CORR=real(fftshift(fft(T_CORR,[],2),2));fv=sum(F_CORR(:,m1:m2),2)/lo;
F_CORR=F_CORR-repmat(fv,1,ntraj);%F_CORR(:,end),1,ntraj);
F_CORR=F_CORR./repmat(max(F_CORR,[],2),1,ntraj);% F_CORR_CO=F_CORR;F_CORR_CO(abs(F_CORR_CO)<0.005)=0;
CMass=abs(sum(OMEGA(:,m0:m).*F_CORR(:,m0:m),2)./sum(F_CORR(:,m0:m),2));
figure;plot(srange,CMass,'r-','linewidth',2);axis square; xlabel('Window Start Time (ps)');
ylabel('<\omega>');title('Center of Mass Trajectory');set(gca,'fontsize',labfont,'linewidth',lwidax);
figure;contourf(omega,srange,F_CORR,25,'edgecolor','none');axis square;xlim([0 20]);colormap('jet');colorbar;
xlabel('\omega (cm^{-1})');ylabel('Window Start Time (ps)');title('C(\omega)');set(gca,'fontsize',labfont,'linewidth',lwidax);
figure;contourf(trange,srange,T_CORR,25,'edgecolor','none');axis square;colormap('jet');colorbar;
xlabel('t (ps)');ylabel('Window Start Time (ps)');title('C(t)');set(gca,'fontsize',labfont,'linewidth',lwidax);caxis([-0.5 1.0]);

dist_n=fitdist(CMass,'Normal');
pdf_n=pdf(dist_n,omega);
figure;h1=histogram(CMass,25,'Normalization','pdf','FaceAlpha',1.0,'FaceColor',[1 0 0],'edgecolor',[0 0 0]);
axis square;xlabel('\omega (cm^{-1})');ylabel('P(\omega)');
title('Probability Distribution <\omega>');set(gca,'fontsize',labfont,'linewidth',lwidax);
hold on;plot(omega,pdf_n,'k-','linewidth',4);xlim([0 20]);hold off;
p_n=normpdf(CMass,dist_n.mu,dist_n.sigma);
logL_n=sum(log(p_n));
[nn,~,nb]=histcounts(CMass,40,'Normalization','pdf');
AIC_n=2*2-2*logL_n;
BIC_n=-2*logL_n+2*log(200);
chi2_n=sum((nn(nb)'-p_n).^2./p_n);

dist_ln=fitdist(CMass,'Lognormal');
pdf_ln=pdf(dist_ln,omega);
figure;h2=histogram(CMass,25,'Normalization','pdf','FaceAlpha',1.0,'FaceColor',[1 0 0],'edgecolor',[0 0 0]);
axis square;xlabel('\omega (cm^{-1})');ylabel('P(\omega)');
title('Probability Distribution <\omega>');set(gca,'fontsize',labfont,'linewidth',lwidax);
hold on;plot(omega,pdf_ln,'k-','linewidth',4);xlim([0 20]);hold off;
p_ln=lognpdf(CMass,dist_ln.mu,dist_ln.sigma);
logL_ln=sum(log(p_ln));
[lnn,~,lnb]=histcounts(CMass,40,'Normalization','pdf');
AIC_ln=2*2-2*logL_ln;
BIC_ln=-2*logL_ln+2*log(200);
chi2_ln=sum((lnn(lnb)'-p_ln).^2./p_ln);

options = statset('MaxIter',10000, 'MaxFunEvals',10000);
GMM=fitgmdist(CMass,2,'Options',options,'CovarianceType','full');
pdf_normmixture = @(x,p,mu1,mu2,sigma1,sigma2) ...
                         p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);                 
paramEsts=[GMM.ComponentProportion(1) GMM.mu(1) GMM.mu(2) sqrt(GMM.Sigma(1)) sqrt(GMM.Sigma(2))];
pdfgrid = pdf_normmixture(omega,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
figure;h3=histogram(CMass,25,'Normalization','pdf','FaceAlpha',1.0,'FaceColor',[1 0 0],'edgecolor',[0 0 0]);
axis square;xlabel('\omega (cm^{-1})');ylabel('P(\omega)');
title('Probability Distribution <\omega>');set(gca,'fontsize',labfont,'linewidth',lwidax);
hold on;plot(omega,pdfgrid,'k-','linewidth',4);xlim([0 20]);hold off;
[Gn,~,Gb]=histcounts(CMass,40,'Normalization','pdf');
p_G=pdf_normmixture(CMass,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
AIC_G=GMM.AIC;
BIC_G=GMM.BIC;
chi2_G=sum((Gn(Gb)'-p_G).^2./p_G);

AIC_all=[AIC_n AIC_ln AIC_G];
BIC_all=[BIC_n BIC_ln BIC_G];
chi2_all=[chi2_n chi2_ln chi2_G];

% dlmwrite('CMass_mbd.txt',CMass);dlmwrite('params_mbd.txt',paramEsts);dlmwrite('omega.txt',omega);
% dlmwrite('CMass_mbc.txt',CMass);dlmwrite('params_mbc.txt',paramEsts);dlmwrite('omega.txt',omega);
% figure;mesh(trange,srange,T_CORR);axis square;colormap('jet');colorbar;
% figure;histogram(K);figure;histogram(TAU,1000);


% pdf_normmixture = @(x,p,mu1,mu2,sigma1,sigma2) ...
%                          p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);                 
% pStart=.5;muStart=[50 150];sigmaStart=[25 25];start=[pStart muStart sigmaStart];
% lb=[0 0 0 0 0];ub=[1 300 300 100 100];options=statset('MaxIter',1000, 'MaxFunEvals',2000);
% paramEsts = mle(CMass, 'pdf',pdf_normmixture, 'start',start,'lower',lb, 'upper',ub, 'options',options);
% pdfgrid = pdf_normmixture(omega,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
% figure;h3=histogram(CMass,50,'Normalization','pdf','FaceAlpha',1.0,'FaceColor',[1 0 0],'edgecolor',[0 0 0]);
% axis square;xlabel('\omega (cm^{-1})');ylabel('P(\omega)');
% title('Probability Distribution <\omega>');set(gca,'fontsize',labfont,'linewidth',lwidax);
% hold on;plot(omega,pdfgrid,'k-','linewidth',4);xlim([0 300]);hold off;
