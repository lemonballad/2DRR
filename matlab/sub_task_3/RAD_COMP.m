hr90=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an90\E_SS\hRAD_90.txt');
hr95=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an95\E_SS\hRAD_95.txt');
hr100=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an100\E_SS\hRAD_100.txt');
hr105=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an105\E_SS\hRAD_105.txt');
hr110=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an110\E_SS\hRAD_110.txt');
hr120=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an120\E_SS\hRAD_120.txt');

pr90=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an90\E_SS\pRAD_90.txt');
pr95=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an95\E_SS\pRAD_95.txt');
pr100=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an100\E_SS\pRAD_100.txt');
pr105=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an105\E_SS\pRAD_105.txt');
pr110=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an110\E_SS\pRAD_110.txt');
pr120=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an120\E_SS\pRAD_120.txt');

hr=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an120\E_SS\R_h_120.txt');
pr=dlmread('C:\Users\Thomas\Desktop\MBC_FE\an120\E_SS\R_p_120.txt');

labfont=16;lwidax=2;

% figure;plot(hr,hr90,hr,hr95,hr,hr100,hr,hr105,hr,hr110,hr,hr120);
% set(gca,'fontsize',labfont,'linewidth',lwidax);
% figure;plot(pr,pr90,pr,pr95,pr,pr100,pr,pr105,pr,pr110,pr,pr120);
% set(gca,'fontsize',labfont,'linewidth',lwidax);

x0=250;
figure;plot(hr(x0:end),hr90(x0:end),hr(x0:end),hr95(x0:end),...
    hr(x0:end),hr100(x0:end),hr(x0:end),hr105(x0:end),...
    hr(x0:end),hr110(x0:end),hr(x0:end),hr120(x0:end));
set(gca,'fontsize',labfont,'linewidth',lwidax);

x0=100;
figure;plot(pr(x0:end),pr90(x0:end),pr(x0:end),pr95(x0:end),...
    pr(x0:end),pr100(x0:end),pr(x0:end),pr105(x0:end),...
    pr(x0:end),pr110(x0:end),pr(x0:end),pr120(x0:end));
set(gca,'fontsize',labfont,'linewidth',lwidax);

