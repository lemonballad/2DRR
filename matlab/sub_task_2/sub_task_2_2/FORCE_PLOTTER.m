clear all;
%% Read in data from files.
prot_40s=dlmread('C:\Users\Thomas\Desktop\FORCE_PLOTTER\prot_40s.txt',',');
prot_60s=dlmread('C:\Users\Thomas\Desktop\FORCE_PLOTTER\prot_60s.txt',',');
prot_90s=dlmread('C:\Users\Thomas\Desktop\FORCE_PLOTTER\prot_90s.txt',',');
delta_r=dlmread('C:\Users\Thomas\Desktop\FORCE_PLOTTER\delta_r.txt',',');
RGB_COMb=Tiff('C:\Users\Thomas\Desktop\FORCE_PLOTTER\RGB_COMb.tif');
[RGB_COMb,a]=RGB_COMb.readRGBAImage();

%% Calculate magnitude of forces on groups of residues.
f_mag_40s=sqrt(sum(prot_40s.^2,2));
f_mag_60s=sqrt(sum(prot_60s.^2,2));
f_mag_90s=sqrt(sum(prot_90s.^2,2));

%% Prepare points for image corners.
mag_r=delta_r(end)-delta_r(1);
mag_f=2500;
[iy,ix]=size(RGB_COMb);
y0=2950;yf=2050;dy=(y0-yf)/(iy-1);dx=iy*dy/ix*mag_r/mag_f;
x0=delta_r(2);xf=dx*(ix-1)+x0;

%% Plot data
labfont=16;lwidax=4;lwid=4;

figure;plot(delta_r,f_mag_40s,'b-',...
    delta_r,f_mag_60s,'r-',...
    delta_r,f_mag_90s,'g-',...
    'linewidth',lwid);
axis square;ylabel('Force (kJ mol^{-1} nm^{-1})');xlabel('\Deltar (nm)');
set(gca,'fontsize',labfont,'linewidth',lwidax);ylim([500 3000]);
hold on;image([x0 xf],[y0 yf],RGB_COMb, 'AlphaData', a);

