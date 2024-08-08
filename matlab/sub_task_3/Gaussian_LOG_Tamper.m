path_source='C:\Users\Thomas\Desktop';
path_dest='C:\Users\Thomas\Desktop';
filename_source='O2Mb_Test';
filename_dest='O2Mb_Test_dest';
file_ext_log='.log';
file_source_log=[path_source '\' filename_source file_ext_log];
file_dest_log=[path_dest '\' filename_dest file_ext_log];

geo=dlmread('C:\Users\Thomas\Desktop\SDs\mbc\mbc_geo_2.txt','');
% geo=geo-repmat(geo(1,:),49,1);

nm=dlmread('C:\Users\Thomas\Desktop\SDs\mbc\mbc_NM_2.txt','');

line_start_geo='                         Standard orientation:';
line_stop_geo=' Leave Link  202 at Thu Dec 15 12:08:30 2016, MaxMem=  2147483648 cpu:         0.1';
preface_geo='(?<preface>\s+\d+\s+\d+\s+\d+\s\s\s\s)';

line_start_nm=' Frequencies --      1.0000';
line_stop_nm=' - Thermochemistry -';
preface_nm='(?<preface>\s+\d+\s+\d+\s\s)';

f_in=fopen(file_source_log,'r+');
f_in_line=fgetl(f_in);
f_out_line=[sprintf('\n')];
read_flag_geo=0;
read_flag_nm=0;
in=0;
ir=0;

while ischar(f_in_line) % Loop over lines in file
    if regexp(f_in_line,line_start_geo) % Read Standard Geometry
        read_flag_geo=1;
        f_out_line=[f_out_line f_in_line sprintf('\n')];  %#ok<*AGROW>
    elseif regexp(f_in_line, line_stop_geo) % Stop reading Standard Geometry
        read_flag_geo=0;
        f_out_line=[f_out_line f_in_line sprintf('\n')];
    elseif regexp(f_in_line,line_start_nm) % Read normal mode
        read_flag_nm=1;
        f_out_line=[f_out_line f_in_line sprintf('\n')];  %#ok<*AGROW>
    elseif regexp(f_in_line, line_stop_nm) % Stop reading normal mode
        read_flag_nm=0;
        f_out_line=[f_out_line f_in_line sprintf('\n')];
    elseif read_flag_geo % Replace geometry data
        pre=regexp(f_in_line,preface_geo,'names');
        if size(pre,1)>0
            ir=ir+1;
            replacement=sprintf('%12.6f%12.6f%12.6f\n',geo(ir,:));
            replacement_line=[pre.preface replacement];
            pre=[];
        else
            replacement_line=[f_in_line sprintf('\n')];
        end
        f_out_line=[f_out_line replacement_line];
    elseif read_flag_nm % Replace normal mode data
        pre=regexp(f_in_line,preface_nm,'names');
        if size(pre,1)>0
            in=in+1;
            replacement=sprintf('%7.2f%7.2f%7.2f\n',nm(ir,:));
            replacement_line=[pre.preface replacement];
            pre=[];
        else
            replacement_line=[f_in_line sprintf('\n')];
        end
        f_out_line=[f_out_line replacement_line];
    else % Non changing lines
        f_out_line=[f_out_line f_in_line sprintf('\n')];
    end % All data has been collected from this file
    f_in_line=fgetl(f_in);
end % loop over lines in file L

f_out=fopen(file_dest_log,'w+');
fprintf(f_out,'%s',f_out_line);
fclose('all');

