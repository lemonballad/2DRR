clear all;
%% Read in data from files.
prot_40s=dlmread('prot_40s.txt',',');
prot_60s=dlmread('prot_60s.txt',',');
prot_90s=dlmread('prot_90s.txt',',');
delta_r=dlmread('delta_r.txt',',');
RGB_COMb=Tiff('RGB_COMb.tif');
[RGB_COMb,a]=RGB_COMb.readRGBAImage();

%% Calculate magnitude of forces on groups of residues.
f_mag_40s=sqrt(sum(prot_40s.^2,2));
f_mag_60s=sqrt(sum(prot_60s.^2,2));
f_mag_90s=sqrt(sum(prot_90s.^2,2));
% Convert kJ/mol/nm to eV/nm
f_mag_40s=f_mag_40s*0.01;
f_mag_60s=f_mag_60s*0.01;
f_mag_90s=f_mag_90s*0.01;

%% Prepare points for image corners.
mag_r=delta_r(end)-delta_r(1);
mag_f=25;
[iy,ix]=size(RGB_COMb);
y0=29.5;yf=20.5;dy=(y0-yf)/(iy-1);dx=iy*dy/ix*mag_r/mag_f;
x0=delta_r(2);xf=dx*(ix-1)+x0;

%% Plot data
labfont=16;lwidax=4;lwid=4;

figure;plot(delta_r,f_mag_40s,'b-',...
    delta_r,f_mag_60s,'r-',...
    delta_r,f_mag_90s,'g-',...
    [-0.024 -0.024],[5 30],'k:',...% Equilibrium distance for COMb
    [-0.036 -0.036],[5 30],'k:',...% Equilibrium distance for deoxy
    'linewidth',lwid);
axis square;ylabel('Force (eV nm^{-1})');xlabel('\Deltar (nm)');
set(gca,'fontsize',labfont,'linewidth',lwidax);
ylim([5.0 30.0]);xlim([delta_r(1) delta_r(end)]);
hold on;image([x0 xf],[y0 yf],RGB_COMb, 'AlphaData', a);

