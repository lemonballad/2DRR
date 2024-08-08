file0='C:\Users\Thomas\Desktop\MBC_FOLDER\MBC_GRO\pre_force.gro';
pattern='\s+(?<a>\d+)(?<b>\w+)\s+\w+\d*\s+(?<c>\d+)';
last_res=0;
fin=fopen(file0);
line=fgetl(fin);
while ischar(line)
    if regexp(line,pattern);
         data=regexp(line,pattern,'names');
         nres=str2double(data.a);
         if last_res~=nres
             last_res=nres;
             res(nres)={data.b};
             ires(nres)=str2double(data.c);
         end
         fres(nres)=str2double(data.c);
         if fres(nres)==1586,break,end;
    end
    line=fgetl(fin);
end
fclose('all');

ires(ires==0)=[];
fres(fres==0)=[];
