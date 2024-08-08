angle=100;
angle_str=['an' num2str(angle)];
nheme_file=1;
nprot_file=41;
nwat_file=41;
nt=50001;
nh=49;
np=1537;
nw=19119;
path=['C:\Users\Thomas\Desktop\MBC_FE\' angle_str '\'];
path_heme=[path 'heme\'];
path_prot=[path 'prot\'];
path_wat=[path 'wat\'];
file_mat_pre='mbc_XYZ';
file_mat_suf='.mat';

%% Read files
%% HEME
full_file_heme_mat=[path_heme file_mat_pre num2str(angle) file_mat_suf];
file_HEME_X=[path_heme 'mbc_X' num2str(angle) '.txt'];
file_HEME_Y=[path_heme 'mbc_Y' num2str(angle) '.txt'];
file_HEME_Z=[path_heme 'mbc_Z' num2str(angle) '.txt'];
filePDB=[path_heme 'heme_' num2str(angle) '_.pdb'];
heme_carts_str='ATOM\s+\d+\s+\S+\s+HEM\s+155\s+(?<x>[-]*\d+[.]+\d+)\s+(?<y>[-]*\d+[.]+\d+)\s+(?<z>[-]*\d+[.]+\d+)\s+';

cHEM_x=zeros(nh,nt,'double');
cHEM_y=zeros(nh,nt,'double');
cHEM_z=zeros(nh,nt,'double');
iatom=0;
itime=1;

mHEME=matfile(full_file_heme_mat,'Writable',true);
mHEME.cHEM_x=[];mHEME.cHEM_y=[];mHEME.cHEM_z=[];
f_in=fopen(filePDB);
line=fgetl(f_in);
HEME_X=fopen(file_HEME_X,'w');
HEME_Y=fopen(file_HEME_Y,'w');
HEME_Z=fopen(file_HEME_Z,'w');

while ischar(line)
    if regexp(line,'ATOM');
        iatom=iatom+1;
        carts_data=regexp(line,heme_carts_str,'names');
        cHEM_x(iatom,itime)=str2double(carts_data.x);
        cHEM_y(iatom,itime)=str2double(carts_data.y);
        cHEM_z(iatom,itime)=str2double(carts_data.z);
        if iatom==nh
            iatom=0;
            itime=itime+1;
        end;
        if mod(itime,100)==0 && iatom==1;
            itime
        end;
    end
    line=fgetl(f_in);
end
fclose(f_in);
mHEME.cHEM_x=cHEM_x;
mHEME.cHEM_y=cHEM_y;
mHEME.cHEM_z=cHEM_z;
data=mHEME.cHEM_x;
fprintf(HEME_X,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
data=mHEME.cHEM_y;
fprintf(HEME_Y,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
data=mHEME.cHEM_z;
fprintf(HEME_Z,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));

%% PROT
full_file_prot_mat=[path_prot file_mat_pre file_mat_suf];
file_PROT_X=[path_prot 'mbc_X' num2str(angle) '.txt'];
file_PROT_Y=[path_prot 'mbc_Y' num2str(angle) '.txt'];
file_PROT_Z=[path_prot 'mbc_Z' num2str(angle) '.txt'];
prot_carts_str='ATOM\s+\d+\s+\S+\s+\S+\s+\d+\s+(?<x>[-]*\d+[.]+\d+)\s+(?<y>[-]*\d+[.]+\d+)\s+(?<z>[-]*\d+[.]+\d+)\s+';

mPROT=matfile(full_file_prot_mat,'Writable',true);
mPROT.cPROT_x=[];mPROT.cPROT_y=[];mPROT.cPROT_z=[];
PROT_X=fopen(file_PROT_X,'w');
PROT_Y=fopen(file_PROT_Y,'w');
PROT_Z=fopen(file_PROT_Z,'w');
parfor ifile=0:nprot_file-1
    M=matfile([path_prot file_mat_pre num2str(ifile) file_mat_suf],'Writable',true);
    M.cPROT_x=[];M.cPROT_y=[];M.cPROT_z=[];
    filePDB=[path_prot 'prot_' num2str(angle) '_' num2str(ifile) '.pdb'];
    f_in=fopen(filePDB);
    line=fgetl(f_in);
    cPROT_x=zeros(np,1500,'double');
    cPROT_y=zeros(np,1500,'double');
    cPROT_z=zeros(np,1500,'double');
    iatom=0;
    itime=1;
    while ischar(line)
        if regexp(line,'ATOM');
            iatom=iatom+1;
            carts_data=regexp(line,prot_carts_str,'names');
            cPROT_x(iatom,itime)=str2double(carts_data.x);
            cPROT_y(iatom,itime)=str2double(carts_data.y);
            cPROT_z(iatom,itime)=str2double(carts_data.z);
            if iatom==np
                iatom=0;
                itime=itime+1;
            end;
            if mod(itime,100)==0 && iatom==1;
                [ifile itime]
            end;
        end
        line=fgetl(f_in);
    end
    fclose(f_in);
    cPROT_x(:,itime:end)=[];
    cPROT_y(:,itime:end)=[];
    cPROT_z(:,itime:end)=[];
    M.cPROT_x=[M.cPROT_x cPROT_x];
    M.cPROT_y=[M.cPROT_y cPROT_y];
    M.cPROT_z=[M.cPROT_z cPROT_z];
end
for ifile=0:nprot_file-1
    M=matfile([path_prot file_mat_pre num2str(ifile) file_mat_suf],'Writable',true);
    mPROT.cPROT_x=[mPROT.cPROT_x M.cPROT_x];
    mPROT.cPROT_y=[mPROT.cPROT_y M.cPROT_y];
    mPROT.cPROT_z=[mPROT.cPROT_z M.cPROT_z];

    ifile
end
data=mPROT.cPROT_x;
fprintf(PROT_X,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
data=mPROT.cPROT_y;
fprintf(PROT_Y,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
data=mPROT.cPROT_z;
fprintf(PROT_Z,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));

%% H2O
full_file_wat_mat=[path_wat file_mat_pre file_mat_suf];
file_H2O_X=[path_wat 'mbc_X' num2str(angle) '.txt'];
file_H2O_Y=[path_wat 'mbc_Y' num2str(angle) '.txt'];
file_H2O_Z=[path_wat 'mbc_Z' num2str(angle) '.txt'];
wat_carts_str='ATOM\s+\d+\s+\S+\s+\S+\s+\d+\s+(?<x>[-]*\d+[.]+\d+)\s+(?<y>[-]*\d+[.]+\d+)\s+(?<z>[-]*\d+[.]+\d+)\s+';

mWAT=matfile(full_file_wat_mat,'Writable',true);
mWAT.cH2O_x=[];mWAT.cH2O_y=[];mWAT.cH2O_z=[];
H2O_X=fopen(file_H2O_X,'w');
H2O_Y=fopen(file_H2O_Y,'w');
H2O_Z=fopen(file_H2O_Z,'w');
for ifile=0:nwat_file-1
    M=matfile([path_wat file_mat_pre num2str(ifile) file_mat_suf],'Writable',true);
    M.cH2O_x=[];M.cH2O_y=[];M.cH2O_z=[];
    filePDB=[path_wat 'wat_' num2str(angle) '_' num2str(ifile) '.pdb'];
    f_in=fopen(filePDB);
    line=fgetl(f_in);
    cH2O_x=zeros(nw,1500,'double');
    cH2O_y=zeros(nw,1500,'double');
    cH2O_z=zeros(nw,1500,'double');
    iatom=0;
    itime=1;
    while ischar(line)
        if regexp(line,'ATOM');
            iatom=iatom+1;
            carts_data=regexp(line,wat_carts_str,'names');
            cH2O_x(iatom,itime)=str2double(carts_data.x);
            cH2O_y(iatom,itime)=str2double(carts_data.y);
            cH2O_z(iatom,itime)=str2double(carts_data.z);
            if iatom==nw
                iatom=0;
                itime=itime+1;
            end;
            if mod(itime,100)==0 && iatom==1;
                [ifile itime]
            end;
        end
        line=fgetl(f_in);
    end
    fclose(f_in);
    cH2O_x(:,itime:end)=[];
    cH2O_y(:,itime:end)=[];
    cH2O_z(:,itime:end)=[];
    M.cH2O_x=[M.cH2O_x cH2O_x];
    M.cH2O_y=[M.cH2O_y cH2O_y];
    M.cH2O_z=[M.cH2O_z cH2O_z];
end
for ifile=0:nwat_file-1
    M=matfile([path_wat file_mat_pre num2str(ifile) file_mat_suf],'Writable',true);
    mWAT.cH2O_x=[mWAT.cH2O_x M.cH2O_x];
    mWAT.cH2O_y=[mWAT.cH2O_y M.cH2O_y];
    mWAT.cH2O_z=[mWAT.cH2O_z M.cH2O_z];
    ifile
end
    data=mWAT.cH2O_x;
    fprintf(H2O_X,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    data=mWAT.cH2O_y;
    fprintf(H2O_Y,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    data=mWAT.cH2O_z;
    fprintf(H2O_Z,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));

