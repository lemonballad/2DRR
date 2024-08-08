file0='C:\Users\Thomas\Desktop\MBC_GAU_GJF\COHEME.gjf';
str_template='(?<starter>\s+\w+\s+)(?<x>[-]*\d+[.]+\d+)(?<y>\s+[-]*\d+[.]+\d+)\s+(?<z>[-]*\d+[.]+\d+)';

for ii=1:25
    gjf_out=['C:\Users\Thomas\Desktop\MBC_GAU_GJF\HEME_' num2str(ii) '.gjf'];
    f0=fopen(file0);
    out=fopen(gjf_out,'w+');
    line=fgetl(f0);
    nl=0;
    atoms=[1 48:49 78:96];
    new_text=[];
    while ischar(line)
        if regexp(line,str_template)
            nl=nl+1;
            data=regexp(line,str_template,'names');
            vec0=[str2double(data.x) str2double(data.y) str2double(data.z)];
            if ~isempty(find(atoms==nl,1))
                vec=vec0+[0 0 (ii-16)*0.1];
                new_line=[data.starter data.x data.y sprintf('%14.8f',vec(3))];
                new_text=[new_text sprintf('%s\n',new_line)];
            else
                new_text=[new_text sprintf('%s\n',line)];
            end
        else
            new_text=[new_text sprintf('%s\n',line)];
        end
        line=fgetl(f0);
    end
    new_text=regexprep(new_text,'%chk=COHEME.chk',['%chk=HEME_' int2str(ii) '.chk']);
    fprintf(out,'%s',new_text);
    fclose('all');
end


