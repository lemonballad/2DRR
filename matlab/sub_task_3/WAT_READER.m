angle=90;
angle_str=['an' num2str(angle)];
nwat_file=41;
nt=50001;
nw=19119;
path=['C:\Users\Thomas\Desktop\MBC_FE\' angle_str '\'];
path_wat=[path 'wat\'];
file_mat_pre='wat_XYZ';
file_mat_suf='.mat';

file_H2O_X=[path_wat 'wat_X' num2str(angle) '.txt'];
file_H2O_Y=[path_wat 'wat_Y' num2str(angle) '.txt'];
file_H2O_Z=[path_wat 'wat_Z' num2str(angle) '.txt'];

H2O_X=fopen(file_H2O_X,'w');
H2O_Y=fopen(file_H2O_Y,'w');
H2O_Z=fopen(file_H2O_Z,'w');

for ifile=0:nwat_file-1
    mWAT=matfile([path_wat file_mat_pre num2str(ifile) file_mat_suf],'Writable',true);
    data=mWAT.cH2O_x;
    fprintf(H2O_X,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    data=mWAT.cH2O_y;
    fprintf(H2O_Y,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    data=mWAT.cH2O_z;
    fprintf(H2O_Z,'%8.2f%8.2f%8.2f%8.2f%8.2f\n',data(:));
    ifile
end
