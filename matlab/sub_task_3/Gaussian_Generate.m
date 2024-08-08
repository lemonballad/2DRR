path_source='C:\Users\Thomas\Desktop\GaussView 3.09\Myoglobin 2016\MBN_F_Input\2-29-16';
path_dest='C:\Users\Thomas\Desktop\GaussView 3.09\Myoglobin 2016\MBN_F_Input\2-29-16';
filename_source='MBN';
filename_prefix_dest_L='MBN_L_';
filename_prefix_dest_R='MBN_R_';
file_ext_gjf='.gjf';
file_ext_chk='.chk';
file_source_gjf=[path_source '\' filename_source file_ext_gjf];
file_source_chk=[path_source '\' filename_source file_ext_chk];
line_old_L='   D37          -81.29974504';
line_old_R='   D13          -81.13323975';
line_new_prefix_L='   D37';
line_new_prefix_R='   D13';
line_new_suffix_L='.00000000';
line_new_suffix_R='.00000000';

file_contents=fileread(file_source_gjf);

for suffix=-90:5:90
    switch 4-length(int2str(suffix))
        case 0; spaces_str='         ';
        case 1; spaces_str='          ';
        case 2; spaces_str='           ';
        case 3; spaces_str='            ';
    end
    
    line_new_L=[line_new_prefix_L spaces_str int2str(-81+suffix) line_new_suffix_L];
    line_new_R=[line_new_prefix_R spaces_str int2str(-81+suffix) line_new_suffix_R];

    file_contents_current_L=regexprep(file_contents,'%chk=MBN.chk',['%chk=MBN_L_' int2str(suffix) '.chk']);
    file_contents_current_R=regexprep(file_contents,'%chk=MBN.chk',['%chk=MBN_R_' int2str(suffix) '.chk']);
    file_contents_current_L=regexprep(file_contents_current_L,line_old_L,line_new_L);
    file_contents_current_R=regexprep(file_contents_current_R,line_old_R,line_new_R);

    file_dest_L_gjf=...
        [path_dest '\' filename_prefix_dest_L int2str(suffix) file_ext_gjf];
    file_dest_R_gjf=...
        [path_dest '\' filename_prefix_dest_R int2str(suffix) file_ext_gjf];
    
    file_dest_L_chk=...
        [path_dest '\' filename_prefix_dest_L int2str(suffix) file_ext_chk];
    file_dest_R_chk=...
        [path_dest '\' filename_prefix_dest_R int2str(suffix) file_ext_chk];
    
    fout_L=fopen(file_dest_L_gjf,'w+');
    fout_R=fopen(file_dest_R_gjf,'w+');

    fprintf(fout_L,'%s',file_contents_current_L);
    fprintf(fout_R,'%s',file_contents_current_R);

    fclose(fout_L);
    fclose(fout_R);
    
    copyfile(file_source_chk,file_dest_L_chk);
    copyfile(file_source_chk,file_dest_R_chk);
end

fclose('all');