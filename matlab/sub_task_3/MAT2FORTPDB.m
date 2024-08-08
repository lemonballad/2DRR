file_HEME_X='C:\Users\Thomas\Desktop\mbd\HEME\mbd_X.txt';
file_HEME_Y='C:\Users\Thomas\Desktop\mbd\HEME\mbd_Y.txt';
file_HEME_Z='C:\Users\Thomas\Desktop\mbd\HEME\mbd_Z.txt';
file_HEME_q='C:\Users\Thomas\Desktop\mbd\HEME\mbd_q.txt';

file_PROT_X='C:\Users\Thomas\Desktop\mbd\PROT\mbd_X.txt';
file_PROT_Y='C:\Users\Thomas\Desktop\mbd\PROT\mbd_Y.txt';
file_PROT_Z='C:\Users\Thomas\Desktop\mbd\PROT\mbd_Z.txt';
file_PROT_q='C:\Users\Thomas\Desktop\mbd\PROT\mbd_q.txt';

file_H2O_X='C:\Users\Thomas\Desktop\mbd\H2O\mbd_X.txt';
file_H2O_Y='C:\Users\Thomas\Desktop\mbd\H2O\mbd_Y.txt';
file_H2O_Z='C:\Users\Thomas\Desktop\mbd\H2O\mbd_Z.txt';
file_H2O_q='C:\Users\Thomas\Desktop\mbc\H2O\mbc_q.txt';

HEME_X=fopen(file_HEME_X,'w');
HEME_Y=fopen(file_HEME_Y,'w');
HEME_Z=fopen(file_HEME_Z,'w');

PROT_X=fopen(file_PROT_X,'w');
PROT_Y=fopen(file_PROT_Y,'w');
PROT_Z=fopen(file_PROT_Z,'w');

H2O_X=fopen(file_H2O_X,'w');
H2O_Y=fopen(file_H2O_Y,'w');
H2O_Z=fopen(file_H2O_Z,'w');

for ifile=0:100
    tic
    mHEME=matfile(['C:\Users\Thomas\Desktop\mbd\HEME\mbd_XYZ' num2str(ifile) '.mat']);
    mPROT=matfile(['C:\Users\Thomas\Desktop\mbd\PROT\mbd_XYZ' num2str(ifile) '.mat']);
    mH2O=matfile(['C:\Users\Thomas\Desktop\mbd\H2O\mbd_XYZ' num2str(ifile) '.mat']);

    data=mHEME.cHEM_x;
    fprintf(HEME_X,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    data=mHEME.cHEM_y;
    fprintf(HEME_Y,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    data=mHEME.cHEM_z;
    fprintf(HEME_Z,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    
    data=mPROT.cPROT_x;
    fprintf(PROT_X,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    data=mPROT.cPROT_y;
    fprintf(PROT_Y,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    data=mPROT.cPROT_z;
    fprintf(PROT_Z,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    
    data=mH2O.cH2O_x;
    fprintf(H2O_X,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    data=mH2O.cH2O_y;
    fprintf(H2O_Y,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    data=mH2O.cH2O_z;
    fprintf(H2O_Z,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    toc
    ifile
end

fclose('all');
% dlmwrite(file_HEME_q,mHEME.q(:,1),'precision','%.2f');
% dlmwrite(file_PROT_q,mPROT.q(:,1),'precision','%.2f');
% dlmwrite(file_H2O_q,mH2O.cH2O_q(:,1),'precision','%.2f');
