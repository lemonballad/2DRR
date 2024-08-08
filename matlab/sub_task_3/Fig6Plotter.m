
t=dlmread('time.txt');omega=dlmread('omega.txt');
mbdSS=dlmread('mbdSS.txt');f_corr_D=dlmread('f_corr_D.txt');
CMass_D=dlmread('CMass_D.txt');pdf_D=dlmread('pdf_D.txt');
mbcSS=dlmread('mbcSS.txt');f_corr_C=dlmread('f_corr_C.txt');
CMass_C=dlmread('CMass_C.txt');pdf_C=dlmread('pdf_C.txt');

labfont=16;lwidax=2;
figure;
subplot(2,3,1);plot(t,mbdSS,'Color',[0 0 1],'Linewidth',lwidax);axis square;set(gca,'fontsize',labfont,'linewidth',lwidax);
xlabel('t (ps)');ylabel('E_{Solute-Solvent} (cm^{-1})');
subplot(2,3,2);semilogy(omega/2/pi,(f_corr_D),'Color',[0 0 1],'Linewidth',lwidax);axis square;set(gca,'fontsize',labfont,'linewidth',lwidax);
xlim([0 3000]);xlabel('\omega (cm^{-1})');ylabel('C(\omega)');
subplot(2,3,3);histogram(CMass_D,20,'Normalization','pdf','FaceAlpha',1.0,'FaceColor',[0 0 1],'edgecolor',[0 0 0]);
axis square;xlabel('\omega (cm^{-1})');ylabel('P(\omega)');set(gca,'fontsize',labfont,'linewidth',lwidax);
hold on;plot(omega,pdf_D,'k-','linewidth',4);xlim([0 20]);hold off;axis square;ylim([0 0.5]);
subplot(2,3,4);plot(t,mbcSS,'Color',[1 0 0],'Linewidth',lwidax);axis square;set(gca,'fontsize',labfont,'linewidth',lwidax);
xlabel('t (ps)');ylabel('E_{Solute-Solvent} (cm^{-1})');ylim([-38000 -28000]);
subplot(2,3,5);semilogy(omega/2/pi,(f_corr_C),'Color',[1 0 0],'Linewidth',lwidax);axis square;set(gca,'fontsize',labfont,'linewidth',lwidax);
xlim([0 3000]);xlabel('\omega (cm^{-1})');ylabel('C(\omega)');
subplot(2,3,6);histogram(CMass_C,20,'Normalization','pdf','FaceAlpha',1.0,'FaceColor',[1 0 0],'edgecolor',[0 0 0]);
axis square;xlabel('\omega (cm^{-1})');ylabel('P(\omega)');set(gca,'fontsize',labfont,'linewidth',lwidax);
hold on;plot(omega,pdf_C,'k-','linewidth',4);xlim([0 20]);hold off;axis square;ylim([0 0.5]);