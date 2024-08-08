angle_str=['an' num2str(angle)];

mHEME=matfile(['C:\Users\Thomas\Desktop\MBC_FE\'...
    angle_str '\heme\heme_XYZ' num2str(angle) '.mat'],'Writable',false);
file_HEME_X=['C:\Users\Thomas\Desktop\MBC_FE\'...
    angle_str '\heme\heme_X' num2str(angle) '.txt'];
file_HEME_Y=['C:\Users\Thomas\Desktop\MBC_FE\'...
    angle_str '\heme\heme_Y' num2str(angle) '.txt'];
file_HEME_Z=['C:\Users\Thomas\Desktop\MBC_FE\'...
    angle_str '\heme\heme_Z' num2str(angle) '.txt'];
HEME_X=fopen(file_HEME_X,'w');
HEME_Y=fopen(file_HEME_Y,'w');
HEME_Z=fopen(file_HEME_Z,'w');
data=mHEME.cHEM_x;
fprintf(HEME_X,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
data=mHEME.cHEM_y;
fprintf(HEME_Y,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
data=mHEME.cHEM_z;
fprintf(HEME_Z,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));

file_PROT_X=['C:\Users\Thomas\Desktop\MBC_FE\'...
    angle_str '\prot\prot_X' num2str(angle) '.txt'];
file_PROT_Y=['C:\Users\Thomas\Desktop\MBC_FE\'...
    angle_str '\prot\prot_Y' num2str(angle) '.txt'];
file_PROT_Z=['C:\Users\Thomas\Desktop\MBC_FE\'...
    angle_str '\prot\prot_Z' num2str(angle) '.txt'];
PROT_X=fopen(file_PROT_X,'w');
PROT_Y=fopen(file_PROT_Y,'w');
PROT_Z=fopen(file_PROT_Z,'w');
for ifile=0:160
    mPROT=matfile(['C:\Users\Thomas\Desktop\MBC_FE\'...
        angle_str '\prot\prot_XYZ' num2str(ifile) '.mat'],'Writable',false);
    data=mPROT.cPROT_x;
    fprintf(PROT_X,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    data=mPROT.cPROT_y;
    fprintf(PROT_Y,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    data=mPROT.cPROT_z;
    fprintf(PROT_Z,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    ifile
end

for angle=[95 105]
file_H2O_X=['C:\Users\Thomas\Desktop\MBC_FE\'...
    angle_str '\wat\' 'wat_X' num2str(angle) '.txt'];
file_H2O_Y=['C:\Users\Thomas\Desktop\MBC_FE\'...
    angle_str '\wat\' 'wat_Y' num2str(angle) '.txt'];
file_H2O_Z=['C:\Users\Thomas\Desktop\MBC_FE\'...
    angle_str '\wat\' 'wat_Z' num2str(angle) '.txt'];

H2O_X=fopen(file_H2O_X,'w');
H2O_Y=fopen(file_H2O_Y,'w');
H2O_Z=fopen(file_H2O_Z,'w');

for ifile=0:160
    mWAT=matfile(['C:\Users\Thomas\Desktop\MBC_FE\'...
        angle_str '\wat\' 'wat_XYZ' num2str(ifile) '.mat'],'Writable',false);
    data=mWAT.cH2O_x;
    fprintf(H2O_X,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    data=mWAT.cH2O_y;
    fprintf(H2O_Y,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    data=mWAT.cH2O_z;
    fprintf(H2O_Z,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    ifile
end
end
