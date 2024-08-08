figure;plot(rd_angle_L_MBN(9:34),freq_L_MBN_mat(9:34,46));
figure;plot(rd_angle_L_MBO(12:34),freq_L_MBO_mat(12:34,46));

file_name_MBN1='C:\users\thomas\desktop\MBN_R_Norm.txt';
file_name_MBO1='C:\users\thomas\desktop\MBO_R_Norm.txt';
file_name_MBN2='C:\users\thomas\desktop\MBN_R_Norm_Adj.txt';
file_name_MBO2='C:\users\thomas\desktop\MBO_R_Norm_Adj.txt';
clear MBN_Angle_Full MBN_Angle MBN_Freq_Full MBN_Freq;
clear MBO_Angle_Full MBO_Angle MBO_Freq_Full MBO_Freq;
MBN_Angle_Full=rd_angle_R_MBN(9:32);
MBN_Angle=rd_angle_R_MBN([9:14 17:21 24:32]);
MBN_Freq_Full=freq_R_MBN_mat(9:32,46);
MBN_Freq=[freq_R_MBN_mat(9:14,46);freq_R_MBN_mat(17,46);freq_R_MBN_mat(18,47);freq_R_MBN_mat([19:21 24:32],46)];
MBO_Angle_Full=rd_angle_R_MBO(9:35);
MBO_Angle=rd_angle_R_MBO(9:34);MBO_Angle([])=[];
MBO_Freq_Full=freq_R_MBO_mat(9:35,46);
MBO_Freq=freq_R_MBO_mat(9:34,46);MBO_Freq([])=[];

figure;plot(MBN_Angle,MBN_Freq,MBN_Angle_Full,MBN_Freq_Full,MBO_Angle_Full,MBO_Freq_Full)

DATA1=[];DATA2=[];
DATA1=[MBN_Angle_Full;MBN_Freq_Full']';
DATA2=[MBN_Angle;MBN_Freq']';
dlmwrite(file_name_MBN1,DATA1,'delimiter','\t');
dlmwrite(file_name_MBN2,DATA2,'delimiter','\t');
DATA1=[];DATA2=[];
DATA1=[MBO_Angle_Full;MBO_Freq_Full']';
DATA2=[MBO_Angle;MBO_Freq']';
dlmwrite(file_name_MBO1,DATA1,'delimiter','\t');
dlmwrite(file_name_MBO2,DATA2,'delimiter','\t');
