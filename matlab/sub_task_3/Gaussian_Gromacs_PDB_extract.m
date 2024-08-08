%% Define paths, FOR MBN
path_source='C:\Users\Thomas\Desktop\Gromacs\Output\4-20-16\MBN';
file_prefix='md_0_2_HEME';
file_ext_pdb='.pdb';
filePDB=[path_source '\' file_prefix file_ext_pdb];

%% Define string expressions to read and flags for reading
test='ATOM   1583  O1D HEM B 154      36.330  42.400  24.100  1.00  0.00           O';
dt=0.02;
tmax=1000;
carts_str='ATOM\s+\d+\s+\S+\s+HEM\s+155\s+(?<x>[-]*\d+[.]+\d+)\s+(?<y>[-]*\d+[.]+\d+)\s+(?<z>[-]*\d+[.]+\d+)\s+';
ter_str='ENDMDL';
L1_str='C2D';
L2_str='C3D';
L3_str='CAD';
L4_str='CBD';
R1_str='C3A';
R2_str='C2A';
R3_str='CAA';
R4_str='CBA';

%% Read file
f_in=fopen(filePDB);
L_prop_carts=zeros(tmax/dt+1,4,3,'double');
R_prop_carts=zeros(tmax/dt+1,4,3,'double');
line=fgetl(f_in);
istep=1;
while ischar(line)
    if regexp(line,L1_str);
        carts_data=regexp(line,carts_str,'names');
        L_prop_carts(istep,1,1:3)=[str2double(carts_data.x) str2double(carts_data.y) str2double(carts_data.z)];
    elseif regexp(line,L2_str);
        carts_data=regexp(line,carts_str,'names');
        L_prop_carts(istep,2,1:3)=[str2double(carts_data.x) str2double(carts_data.y) str2double(carts_data.z)];
    elseif regexp(line,L3_str);
        carts_data=regexp(line,carts_str,'names');
        L_prop_carts(istep,3,1:3)=[str2double(carts_data.x) str2double(carts_data.y) str2double(carts_data.z)];
    elseif regexp(line,L4_str);
        carts_data=regexp(line,carts_str,'names');
        L_prop_carts(istep,4,1:3)=[str2double(carts_data.x) str2double(carts_data.y) str2double(carts_data.z)];
    end
    
    if regexp(line,R1_str);
        carts_data=regexp(line,carts_str,'names');
        R_prop_carts(istep,1,1:3)=[str2double(carts_data.x) str2double(carts_data.y) str2double(carts_data.z)];
    elseif regexp(line,R2_str);
        carts_data=regexp(line,carts_str,'names');
        R_prop_carts(istep,2,1:3)=[str2double(carts_data.x) str2double(carts_data.y) str2double(carts_data.z)];
    elseif regexp(line,R3_str);
        carts_data=regexp(line,carts_str,'names');
        R_prop_carts(istep,3,1:3)=[str2double(carts_data.x) str2double(carts_data.y) str2double(carts_data.z)];
    elseif regexp(line,R4_str);
        carts_data=regexp(line,carts_str,'names');
        R_prop_carts(istep,4,1:3)=[str2double(carts_data.x) str2double(carts_data.y) str2double(carts_data.z)];
    end
    
    if regexp(line,ter_str),istep=istep+1;end
    
    line=fgetl(f_in);
    
end

%% Calculate dihedral angles
clear L_v1 L_v2 L_v3 R_v1 R_v2 R_v3 L_c1 L_c2 L_c12 R_c1 R_c2 R_c12 L_d1 L_d2 R_d1 R_d2 
L_v1(:,:)=L_prop_carts(:,2,:)-L_prop_carts(:,1,:);
L_v2(:,:)=L_prop_carts(:,3,:)-L_prop_carts(:,2,:);
L_v3(:,:)=L_prop_carts(:,4,:)-L_prop_carts(:,3,:);
L_v1=L_v1./repmat(sqrt(sum(L_v1.^2,2)),1,3);
L_v2=L_v2./repmat(sqrt(sum(L_v2.^2,2)),1,3);
L_v3=L_v3./repmat(sqrt(sum(L_v3.^2,2)),1,3);

R_v1(:,:)=R_prop_carts(:,2,:)-R_prop_carts(:,1,:);
R_v2(:,:)=R_prop_carts(:,3,:)-R_prop_carts(:,2,:);
R_v3(:,:)=R_prop_carts(:,4,:)-R_prop_carts(:,3,:);
R_v1=R_v1./repmat(sqrt(sum(R_v1.^2,2)),1,3);
R_v2=R_v2./repmat(sqrt(sum(R_v2.^2,2)),1,3);
R_v3=R_v3./repmat(sqrt(sum(R_v3.^2,2)),1,3);

L_c1=cross(L_v1,L_v2,2);
L_c2=cross(L_v2,L_v3,2);
L_c12=cross(L_c1,L_c2,2);
L_c1=L_c1./repmat(sqrt(sum(L_c1.^2,2)),1,3);
L_c2=L_c2./repmat(sqrt(sum(L_c2.^2,2)),1,3);
L_c12=L_c12./repmat(sqrt(sum(L_c12.^2,2)),1,3);

R_c1=cross(R_v1,R_v2,2);
R_c2=cross(R_v2,R_v3,2);
R_c12=cross(R_c1,R_c2,2);
R_c1=R_c1./repmat(sqrt(sum(R_c1.^2,2)),1,3);
R_c2=R_c2./repmat(sqrt(sum(R_c2.^2,2)),1,3);
R_c12=R_c12./repmat(sqrt(sum(R_c12.^2,2)),1,3);

L_d1=dot(L_c12,L_v2,2);
L_d2=dot(L_c1,L_c2,2);

R_d1=dot(R_c12,R_v2,2);
R_d2=dot(R_c1,R_c2,2);

L_dihedrals=round(atan2(L_d1,L_d2)*180/pi)+81;
R_dihedrals=round(atan2(R_d1,R_d2)*180/pi)+81;

L_dihedrals(L_dihedrals>=180)=L_dihedrals(L_dihedrals>=180)-360;
L_dihedrals(L_dihedrals<-180)=360+L_dihedrals(L_dihedrals<-180);
R_dihedrals(R_dihedrals>=180)=R_dihedrals(R_dihedrals>=180)-360;
R_dihedrals(R_dihedrals<-180)=360+R_dihedrals(R_dihedrals<-180);

%%
freq_PDB=[];
for ii=1:length(L_dihedrals)
    freq_PDB(ii)=freq_map(L_dihedrals(ii)+181,R_dihedrals(ii)+181);
end
frange=1:length(freq_PDB);
figure;plot(L_dihedrals);
figure;plot(R_dihedrals);
figure;plot(freq_PDB(frange));
mean_freq=mean(freq_PDB(frange));
var_freq=var(freq_PDB(frange));
testcorr=freq_PDB(frange);
tempcorr=xcorr(testcorr-mean_freq,'coeff');
% tempcorr=xcorr(testcorr);
testcorr=tempcorr(length(testcorr):end);
time=0:dt:(length(testcorr)-1)*dt;
f=fit(time',testcorr','exp2','StartPoint',[0.6 -1/0.1 0.4 -1/200]);%'Lower',[max(testcorr)-0.001 -50],'Upper',[max(testcorr) 1/max(time)]);
% f=fit(time',testcorr,'exp1');
f1=f.a*exp(f.b*time)+f.c*exp(f.d*time);
f2=exp(f.d*time);
figure;plot(time,testcorr,time,f1,time,testcorr-f1,time,md1);
figure;loglog(time,testcorr);

% modelfun = @(b,x)b(1)*cos(2*pi*x/b(2)+b(3))+b(4)*sin(2*pi*x/b(5)+b(6))+b(7);
% beta0=[0.0004 70 pi/6 -0.005 200 -pi mean(filcorr)];
% md1=-0.001*cos(2*pi*time/47+-3)+-0.003*sin(2*pi*time/175+-6)-0.001;
% modelfun = @(b,x)b(1)*cos(2*pi*x/b(2)+b(3))+b(4);%;
% beta0=[0.04 70 -pi mean(filcorr)];
% mdl = fitnlm(time,filter(0.01*ones(1,100),1,testcorr-f1),modelfun,beta0);
% md1=0.008*cos(2*pi*time/76+-0.18)-0.001;figure;plot(time,testcorr-f1-md1);

% filcorr=filter(0.01*ones(1,100),1,testcorr);figure;plot(time,filcorr);
% ffff=apodfun(time,0,0,100,500,false);
% figure;plot((-1/2:1/length(freq_PDB):1/2-1/length(freq_PDB))/0.02/0.00003,abs(fftshift(fft(ffff.*filcorr))));xlim([-1000 1000]);
% 
% figure;plot((-1/2:1/length(freq_PDB):1/2-1/length(freq_PDB))/0.02/0.00003,abs(fftshift(fft(testcorr-f1))));xlim([-5000 5000]);
% hold on
% plot((-1/2:1/length(freq_PDB):1/2-1/length(freq_PDB))/0.02/0.00003,abs(fftshift(fft(f1))));xlim([-5000 5000]);
% hold on
% plot((-1/2:1/length(freq_PDB):1/2-1/length(freq_PDB))/0.02/0.00003,abs(fftshift(fft(testcorr))));xlim([-5000 5000]);


ffff=apodfun(time,0,0,18,16,false);figure;plot(time,testcorr,time,testcorr.*ffff,time,f1);ylim([-0.04 0.05]);xlim([0 300]);
ff=((-1/2:1/length(freq_PDB):1/2-1/length(freq_PDB))+1/length(freq_PDB)/2)*2*pi/20/0.00003;
fcorr=real(fftshift(fft(ffff.*testcorr)));
ff1=abs(fftshift(fft(f1)));
ff2=real(fftshift(fft(f2)));
f3=real(1./(abs(f.b)+2*pi*1i*ff));
figure;plot(ff,fcorr/max(fcorr),ff,ff2/max(ff2),'Linewidth',3);%(-1/2:1/length(freq_PDB):1/2-1/length(freq_PDB))/0.02/0.00003,abs(fftshift(fft(testcorr))),...
legend('filtered data','fit exponentials');xlim([0 10]);xlabel('\omega (cm^-^1)');

% figure;plot((-1/2:1/length(freq_PDB):1/2-1/length(freq_PDB))/0.02/0.00003,ffff.*testcorr,...
%     (-1/2:1/length(freq_PDB):1/2-1/length(freq_PDB))/0.02/0.00003,testcorr);