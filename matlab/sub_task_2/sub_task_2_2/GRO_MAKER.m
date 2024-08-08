file0='C:\Users\Thomas\Desktop\MBC_GRO\pre_force.gro';
original_txt=fileread(file0);
% old_line1='  155HEM     FE 1538   3.043   3.705   2.677 -0.0426  0.0061 -0.0030';
% old_line2='  155HEM    C1O 1585   3.202   3.632   2.709 -0.3271 -0.7221 -0.2430';
% old_line3='  155HEM    O1C 1586   3.313   3.632   2.697 -0.2775 -0.4848  0.2343';
old_line1='   3.043   3.705   2.677';
old_line2='   3.202   3.632   2.709';
old_line3='   3.313   3.632   2.697';
r0=[3.094+3.000+3.063+3.150...
    3.832+3.568+3.557+3.820...
    2.824+2.816+2.545+2.553]/4;
r1=[3.043 3.705 2.677];
r2=[3.202 3.632 2.709];
r3=[3.313 3.632 2.697];
v1=[3.094-3.000 3.832-3.568 2.824-2.816];
v2=[3.094-3.000 3.832-3.557 2.824-2.545];
n=cross(v1,v2);n=n/norm(n);
dr=0.01;
r=-0.1:dr:0.1;r=r-0.024;
lr=length(r);
for ii=1:lr
    R1=r1+r(ii)*n;
    R2=r2+r(ii)*n;
    R3=r3+r(ii)*n;
    gro_out=['C:\Users\Thomas\Desktop\MBC_GRO\pre_force_' num2str(ii) '.gro'];
    new_line1=[sprintf('%8.3f',R1(1)) sprintf('%8.3f',R1(2)) sprintf('%8.3f',R1(3))];
    new_line2=[sprintf('%8.3f',R2(1)) sprintf('%8.3f',R2(2)) sprintf('%8.3f',R2(3))];
    new_line3=[sprintf('%8.3f',R3(1)) sprintf('%8.3f',R3(2)) sprintf('%8.3f',R3(3))];
    new_txt=original_txt; %#ok<*NASGU>
    new_txt=regexprep(new_txt,old_line1,new_line1);
    new_txt=regexprep(new_txt,old_line2,new_line2);
    new_txt=regexprep(new_txt,old_line3,new_line3);
    out=fopen(gro_out,'w+');
    fprintf(out,'%s',new_txt);
    fclose('all');
end


