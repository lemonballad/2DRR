% t=1:50001;nan_axis=[nan_axis 50001];
mHEME=matfile('C:\Users\Thomas\Desktop\mbc\HEME\mbc_XYZ.mat');
mH2O=matfile('C:\Users\Thomas\Desktop\mbc\H2O\mbc_XYZ.mat');
mPROT=matfile('C:\Users\Thomas\Desktop\mbc\PROT\mbc_XYZ.mat');
cHEM_x=mHEME.cHEM_x(:,1:74970);cHEM_y=mHEME.cHEM_y(:,1:74970);cHEM_z=mHEME.cHEM_z(:,1:74970);
nt=size(cHEM_x,2);
t=0:0.004:(nt-1)*0.004;

% q_temp=M.cH2O_x;
% q=NaN(19122,50001);q(:,t_axis)=q_temp;clear q_temp;q(:,nan_axis)=(q(:,[nan_axis(1:20)-1 49999])+q(:,[nan_axis(1:19)+1 49999 49999]))/2;
q=mH2O.cH2O_x;
Q=zeros(size(q),'double');
for iatom=1:size(cHEM_x,1)
    Q=Q+(q-repmat(cHEM_x(iatom,:),size(q,1),1)).^2;
    [1 iatom]
end
q=mH2O.cH2O_y;
% q_temp=M.cH2O_y;
% q=NaN(19122,50001);q(:,t_axis)=q_temp;clear q_temp;q(:,nan_axis)=(q(:,[nan_axis(1:20)-1 49999])+q(:,[nan_axis(1:19)+1 49999 49999]))/2;
for iatom=1:size(cHEM_y,1)
    Q=Q+(q-repmat(cHEM_y(iatom,:),size(q,1),1)).^2;
    [2 iatom]
end
q=mH2O.cH2O_z;
% q=NaN(19122,50001);q_temp=M.cH2O_z;q=NaN(19122,50001);
% q(:,t_axis)=q_temp;clear q_temp;q(:,nan_axis)=(q(:,[nan_axis(1:20)-1 49999])+q(:,[nan_axis(1:19)+1 49999 49999]))/2;
for iatom=1:size(cHEM_z,1)
    Q=Q+(q-repmat(cHEM_z(iatom,:),size(q,1),1)).^2;
    [3 iatom]
end
Q=sqrt(Q);
q=mH2O.cH2O_q;
epsilon_0=8.85*10^-12; % In SI units
q=(1/4/pi/epsilon_0)*repmat(cHEM_q(iatom,:),size(q,1),1).*q;
Q=q./Q;
Q=sum(Q,1);
clear q;

Solu_Solv=Q;%zeros(1,nt,'double');
cPROT_x=mPROT.cPROT_x;cPROT_y=mPROT.cPROT_y;cPROT_z=mPROT.cPROT_z;cPROT_q=mPROT.cPROT_q;
for iatom=1:size(cHEM_z,1)
Solu_Solv=Solu_Solv+sum((1/4/pi/epsilon_0)*repmat(cHEM_q(iatom,:),size(cPROT_q,1),1).*cPROT_q./...
        sqrt(...
        (cPROT_x-repmat(cHEM_x(iatom,:),size(cPROT_x,1),1)).^2+...
        (cPROT_y-repmat(cHEM_y(iatom,:),size(cPROT_y,1),1)).^2+...
        (cPROT_z-repmat(cHEM_z(iatom,:),size(cPROT_z,1),1)).^2),1);
    iatom
end
clear cPROT_x cPROT_y cPROT_z cPROT_q cHEM_x cHEM_y cHEM_z cHEM_q
Solu_Solv=Solu_Solv*10^9; % Accounting for the denominator needing to go from nanometers to meters
Solu_Solv=Solu_Solv*(1.6*10^-19)^2; % Converting the charge from e- to coulombs
Solu_Solv=Solu_Solv*5.034*10^22; % Convert from joules to wavenumbers

figure;plot(t,Solu_Solv);
% M.Solu_Solv=Solu_Solv;
