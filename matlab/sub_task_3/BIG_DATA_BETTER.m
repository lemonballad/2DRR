mHEME=matfile('C:\Users\Thomas\Desktop\mbc\HEME\mbc_XYZ.mat');
mPROT=matfile('C:\Users\Thomas\Desktop\mbc\PROT\mbc_XYZ.mat');
mH2O=matfile(['C:\Users\Thomas\Desktop\mbc\mbc_XYZ.mat']);
nt=250001;
dt=0.004;
cHEM_x=mHEME.cHEM_x(:,1:nt);cHEM_y=mHEME.cHEM_y(:,1:nt);cHEM_z=mHEME.cHEM_z(:,1:nt);
t=0:dt:(nt-1)*dt;

cHEM_q=mHEME.q;
    cHEM_x(cHEM_q(:,1)==0,:)=[];
    cHEM_y(cHEM_q(:,1)==0,:)=[];
    cHEM_z(cHEM_q(:,1)==0,:)=[];
    cHEM_q(cHEM_q(:,1)==0,:)=[];
cPROT_q=mPROT.q;
    px=mPROT.cPROT_x;py=mPROT.cPROT_y;pz=mPROT.cPROT_z;
    px(cPROT_q(:,1)==0,:)=[];
    py(cPROT_q(:,1)==0,:)=[];
    pz(cPROT_q(:,1)==0,:)=[];
    cPROT_q(cPROT_q(:,1)==0,:)=[];
cH2O_q=mH2O.cH2O_q(:,1);
Q=zeros(1,nt,'double');
parfor it=1:nt
    ifile=floor((it-1)/2500);
    it_H2O=mod(it,2500);if it_H2O==0,it_H2O=2500;end
    mH2O=matfile(['C:\Users\Thomas\Desktop\mbc\H2O\mbc_XYZ' num2str(ifile) '.mat']);
    hx=mH2O.cH2O_x(:,it_H2O);hy=mH2O.cH2O_y(:,it_H2O);hz=mH2O.cH2O_z(:,it_H2O);
    for iatom=1:size(cHEM_x,1)
        Q(it)=Q(it)+sum(cHEM_q(iatom)*cH2O_q./sqrt(...
            (hx-cHEM_x(iatom,it)).^2+...
            (hy-cHEM_y(iatom,it)).^2+...
            (hz-cHEM_z(iatom,it)).^2));
        Q(it)=Q(it)+sum(cHEM_q(iatom)*cPROT_q./sqrt(...
            (px(:,it)-cHEM_x(iatom,it)).^2+...
            (py(:,it)-cHEM_y(iatom,it)).^2+...
            (pz(:,it)-cHEM_z(iatom,it)).^2));
    end
    it
end
epsilon_0=8.85*10^-12; % In SI units
Q=(1/4/pi/epsilon_0)*Q;

Solu_Solv=Q*10^9; % Accounting for the denominator needing to go from nanometers to meters
Solu_Solv=Solu_Solv*(1.6*10^-19)^2; % Converting the charge from e- to coulombs
Solu_Solv=Solu_Solv*5.034*10^22; % Convert from joules to wavenumbers

figure;plot(t,Solu_Solv);
mHEME.Solu_Solv=Solu_Solv;