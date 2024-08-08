file0='C:\Users\Thomas\Desktop\MBC_FOLDER\PDB\COMbH.pdb';
'ATOM      1  N   VAL A   1      -2.900  17.600  15.500                       N';

'ATOM      1  N   VAL     1      12.000  25.590  37.480  1.00  0.37           N';
%     force=([sqrt(sum(pforce.^2,3)) 0*sqrt(sum(hforce.^2,3))]);
%     force=99*((force-min(min(force)))/max(max(force)));
for ii=1:lr
for jj=1:length(fres)-1;
    force(ii,ires(jj):fres(jj))=sqrt(sum(sum(pforce(ii,ires(jj):fres(jj),:),2).^2,3));
end
% force(ii,1538:end)=sqrt(sum(sum(hforce(ii,:,:),2).^2,3));
force(ii,1538:end)=72+0*sqrt(sum(hforce(ii,:,:).^2,3));
end
force=log10(force);

test_pattern='ATOM\s+(?<natom>\d+)';
pattern='(?<a>\w+\s+\d+\s+\w+\s+\w+\s+\d+)\s+(?<x>[-]*\d+[.]+\d+)\s+(?<y>[-]*\d+[.]+\d+)\s+(?<z>[-]*\d+[.]+\d+)\s+\d+[.]\d+\s+\d+[.]\d+\s+(?<b>\w+)';

r0=[30.94+30.00+30.63+31.50...
    38.32+35.68+35.57+38.20...
    28.24+28.16+25.45+25.53]/4;
r1=[30.43 37.05 26.77];
r2=[32.02 36.32 27.09];
r3=[33.13 36.32 26.97];
v1=[30.94-30.00 38.32-35.68 28.24-28.16];
v2=[30.94-30.00 38.32-35.57 28.24-25.45];
n=cross(v1,v2);n=n/norm(n);
dr=0.1;
r=-1:dr:1;
lr=length(r);

atoms=[1538 1585 1586];

pdb_out=['C:\Users\Thomas\Desktop\MBC_FOLDER\PDB\force_TRAJ.pdb'];
fout=fopen(pdb_out,'w+');
for ii=1:lr
% pdb_out=['C:\Users\Thomas\Desktop\MBC_FOLDER\PDB\force_' num2str(ii) '.pdb'];
% fout=fopen(pdb_out,'w+');
fin=fopen(file0);
new_line=['MODEL' sprintf('%9d',num2str(ii))];
fprintf(fout,'%s\n',new_line);
line=fgetl(fin);
while ischar(line)
    if regexp(line,test_pattern)
        line;
        data=regexp(line,test_pattern,'names');
        natom=str2double(data.natom);
        data=regexp(line,pattern,'names');
        pre_line=data.a;
        x=str2double(data.x);
        y=str2double(data.y);
        z=str2double(data.z);
        post_line=data.b;
        if ~isempty(find(atoms==natom,1))
            new_x=x+r(ii)*n(1);
            new_y=y+r(ii)*n(2);
            new_z=z+r(ii)*n(3);
        else
            new_x=x;
            new_y=y;
            new_z=z;
        end
        new_line=[pre_line sprintf('%12.3f',new_x) sprintf('%8.3f',new_y)...
            sprintf('%8.3f',new_z) sprintf('%6.2f',force(ii,natom))...
            sprintf('%6.2f',force(ii,natom)) sprintf('%12s',post_line)];
        fprintf(fout,'%s\n',new_line);
    else
        line;
        fprintf(fout,'%s\n',line);
    end
    line=fgetl(fin);
end
    fclose(fin);
    ii
end

fclose('all');
