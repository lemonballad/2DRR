%% Define paths, FOR MBN
path_source='C:\Users\Thomas\Desktop\Gromacs\Output\4-29-16\MBN';
file_prefix='md_0_1_HEME';
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
if false
freq_PDB=[];
for ii=1:length(L_dihedrals)
    freq_PDB(ii)=freq_map(L_dihedrals(ii)+181,R_dihedrals(ii)+181);
end
frange=20000:length(freq_PDB);
figure;plot(L_dihedrals(frange));
figure;plot(R_dihedrals(frange));
figure;plot(freq_PDB(frange));
mean_freq=mean(freq_PDB(frange));
var_freq=var(freq_PDB(frange));
testcorr=freq_PDB(frange);
tempcorr=xcorr(testcorr-mean_freq);
% tempcorr=xcorr(testcorr);
testcorr=tempcorr(length(testcorr):end);
time=0:dt:(length(testcorr)-1)*dt;
f=fit(time',testcorr','exp2','StartPoint',[0.6 -1/0.1 0.4 -1/11]);%'Lower',[max(testcorr)-0.001 -50],'Upper',[max(testcorr) 1/max(time)]);
% f=fit(time',testcorr,'exp1');
f1=f.a*exp(f.b*time)+f.c*exp(f.d*time);
f2=exp(f.d*time);
figure;plot(time,testcorr,time,f1,time,testcorr-f1);xlabel('t (ps)');ylabel('<\omega(0)|\omega(t)>');
figure;loglog(time,testcorr);xlabel('log_{10} |t|');ylabel('<log_{10} |\omega(0)|\omega(t)>|');

ffff=apodfun(time,0,0,8,8,false);figure;plot(time,testcorr,time,testcorr.*ffff);xlim([0 300]);ylim([-1 5]*10^5);
Nf=length(freq_PDB);
ff=((-1/2:1/length(freq_PDB):1/2-1/length(freq_PDB))+1/length(freq_PDB)/2)*2*pi/20/0.00003;
fcorr=real(fftshift(fft(ffff.*testcorr,Nf)));
ff1=real(fftshift(fft(f1,Nf)));
ff2=real(fftshift(fft(f2,Nf)));
f3=real(1./(abs(f.d)+1i*ff));
figure;plot(ff,fcorr,ff,ff1/max(ff1)*max(fcorr),'Linewidth',3);%(-1/2:1/length(freq_PDB):1/2-1/length(freq_PDB))/0.02/0.00003,abs(fftshift(fft(testcorr))),...
legend('filtered data','fit exponentials');xlim([0 100]);xlabel('\omega (cm^-^1)');ylabel('<\omega(0)|\omega(t)>');
end
if false
    file_num=5;
    var_freq=var(freq_PDB(frange));
    mbn(file_num).LDA=L_dihedrals';
    mbn(file_num).RDA=R_dihedrals';
    mbn(file_num).freq=freq_PDB;
    mbn(file_num).time=time;
    mbn(file_num).omega=ff;
    mbn(file_num).var=var_freq;
    mbn(file_num).fit=f;
    mbn(file_num).spec_den=fcorr;
    mbn(file_num).corr=testcorr;
    mbn(file_num).apod_param=[0 0 8 8];
    mbn(file_num).frange=frange;
end
