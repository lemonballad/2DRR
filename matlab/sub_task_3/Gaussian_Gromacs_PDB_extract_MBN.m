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

