% Define plotting parameters and read data
labfont=16;lwidax=2;
CMass=dlmread('CMass_mbd.txt');omega=dlmread('omega.txt');paramEsts=dlmread('params_mbd.txt');

% Plot Histogram and Distribution for Normal
dist_n=fitdist(CMass,'Normal');
pdf_n=pdf(dist_n,omega);
figure;h1=histogram(CMass,50,'Normalization','pdf','FaceAlpha',1.0,'FaceColor',[0 0 1],'edgecolor',[0 0 0]);
axis square;xlabel('\omega (cm^{-1})');ylabel('P(\omega)');
title('Probability Distribution <\omega>');set(gca,'fontsize',labfont,'linewidth',lwidax);
hold on;plot(omega,pdf_n,'k-','linewidth',4);xlim([0 300]);hold off;
% Calculate Statistics for Normal;
p_n=normpdf(CMass,dist_n.mu,dist_n.sigma);
logL_n=sum(log(p_n));[nn,~,nb]=histcounts(CMass,50,'Normalization','pdf');
AIC_n=2*2-2*logL_n;BIC_n=-2*logL_n+2*log(200);chi2_n=sum((nn(nb)'-p_n).^2./p_n);

% Plot Histogram and Distribution for LogNormal
dist_ln=fitdist(CMass,'Lognormal');
pdf_ln=pdf(dist_ln,omega);
figure;h2=histogram(CMass,50,'Normalization','pdf','FaceAlpha',1.0,'FaceColor',[0 0 1],'edgecolor',[0 0 0]);
axis square;xlabel('\omega (cm^{-1})');ylabel('P(\omega)');
title('Probability Distribution <\omega>');set(gca,'fontsize',labfont,'linewidth',lwidax);
hold on;plot(omega,pdf_ln,'k-','linewidth',4);xlim([0 300]);hold off;
% Calculate Statistics for LogNormal
p_ln=lognpdf(CMass,dist_ln.mu,dist_ln.sigma);
logL_ln=sum(log(p_ln));[lnn,~,lnb]=histcounts(CMass,50,'Normalization','pdf');
AIC_ln=2*2-2*logL_ln;BIC_ln=-2*logL_ln+2*log(200);chi2_ln=sum((lnn(lnb)'-p_ln).^2./p_ln);

% Plot Histogram and Distribution for Bimodal
pdf_normmixture = @(x,p,mu1,mu2,sigma1,sigma2) ...
                         p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);                 
pdfgrid = pdf_normmixture(omega,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
figure;h3=histogram(CMass,50,'Normalization','pdf','FaceAlpha',1.0,'FaceColor',[0 0 1],'edgecolor',[0 0 0]);
axis square;xlabel('\omega (cm^{-1})');ylabel('P(\omega)');
title('Probability Distribution <\omega>');set(gca,'fontsize',labfont,'linewidth',lwidax);
hold on;plot(omega,pdfgrid,'k-','linewidth',4);xlim([0 300]);hold off;
% Calculate Statistics for Bimodal
[Gn,~,Gb]=histcounts(CMass,50,'Normalization','pdf');
p_G=pdf_normmixture(CMass,paramEsts(1),paramEsts(2),paramEsts(3),paramEsts(4),paramEsts(5));
AIC_G=GMM.AIC;BIC_G=GMM.BIC;chi2_G=sum((Gn(Gb)'-p_G).^2./p_G);

% Collect Statistics
AIC_all=[AIC_n AIC_ln AIC_G];BIC_all=[BIC_n BIC_ln BIC_G];chi2_all=[chi2_n chi2_ln chi2_G];
RNames={'Normal';'LogNor';'Bimodl'};VNames={'AIC';'BIC';'chi2'};

% Print table of statistics
T=table(round(AIC_all*10)'/10,round(BIC_all*10)'/10,round(chi2_all*100)'/100,'RowNames',RNames,'VariableNames',VNames)