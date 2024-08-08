R_flag=0;
L_flag=1;
DATA=[];
file_name='C:\Users\Thomas\Desktop\MBN_L.txt';
% file_name='C:\Users\Thomas\Desktop\MBN_R.txt';
red_mass_stretch=6;
red_mass_bend=2;
u_kg=1.66*10^-27;
c_cms=3*10^10;
mDynPerAng_PerS2=10^2;

dihedral_angle=75+74+37;bond_1=35;bond_2=36;angle=75+35;
% dihedral_angle=75+74+13;bond_1=10;bond_2=11;angle=75+10;

if L_flag
    DATA(:,1)=rd_angle_L_MBN;DATA(:,2)=frc_cnst_L_MBN(dihedral_angle,bond_1,:);
    DATA(:,3)=frc_cnst_L_MBN(dihedral_angle,bond_2,:);DATA(:,4)=frc_cnst_L_MBN(dihedral_angle,angle,:);
end

if R_flag
    DATA(:,1)=rd_angle_R_MBN;DATA(:,2)=frc_cnst_R_MBN(dihedral_angle,bond_1,:);
    DATA(:,3)=frc_cnst_R_MBN(dihedral_angle,bond_2,:);DATA(:,4)=frc_cnst_R_MBN(dihedral_angle,angle,:);
end

DATA(:,2)=1/c_cms*sqrt(DATA(:,2)*mDynPerAng_PerS2/red_mass_stretch/u_kg);
DATA(:,3)=1/c_cms*sqrt(DATA(:,3)*mDynPerAng_PerS2/red_mass_stretch/u_kg);
DATA(:,4)=1/c_cms*sqrt(DATA(:,4)*mDynPerAng_PerS2/red_mass_bend/u_kg);

dlmwrite(file_name,DATA,'delimiter','\t');