path_source='C:\Users\Thomas\Desktop\GaussView 3.09\Myoglobin 2016\MBO_Grid_Input\3-30-16';
path_dest='C:\Users\Thomas\Desktop\GaussView 3.09\Myoglobin 2016\MBO_Grid_Input\3-30-16';
filename_source='MBO';
filename_prefix_dest='MBO_';
% path_source='C:\Users\Thomas\Desktop\GaussView 3.09\Myoglobin 2016\MBN_Grid_Input\3-30-16';
% path_dest='C:\Users\Thomas\Desktop\GaussView 3.09\Myoglobin 2016\MBN_Grid_Input\3-30-16';
% filename_source='MBN';
% filename_prefix_dest='MBN_';
file_ext_gjf='.gjf';
% file_ext_chk='.chk';
file_source_gjf=[path_source '\' filename_source file_ext_gjf];
% file_source_chk=[path_source '\' filename_source file_ext_chk];
line_old_L='   D39          -94.37048385'; %MBO
line_old_R='   D15         -109.09904387'; %MBO
line_new_prefix_L='   D39'; %MBO
line_new_prefix_R='   D15'; %MBO
% line_old_L='   D37          -81.29974504'; %MBN
% line_old_R='   D13          -81.13323975'; %MBN
% line_new_prefix_L='   D37'; %MBN
% line_new_prefix_R='   D13'; %MBN

file_contents=fileread(file_source_gjf);

for suffix_L=-180:5:175              % FOR MBO
    for suffix_R=-180:5:175          % FOR MBO
% for suffix_L=-180:5:175           % FOR MBN
%     for suffix_R=-180:5:175       % FOR MBN
        line_new_L=sprintf('%-14s%14.8f',line_new_prefix_L,(-94.37048385+suffix_L)); %MBO
        line_new_R=sprintf('%-14s%14.8f',line_new_prefix_R,(-109.09904387+suffix_R)); %MBO
%         line_new_L=sprintf('%-14s%14.8f',line_new_prefix_L,(-81.29974504+suffix_L)); %MBN
%         line_new_R=sprintf('%-14s%14.8f',line_new_prefix_R,(-81.13323975+suffix_R)); %MBN
        
        file_contents_current=regexprep(file_contents,'%chk=MBO.chk',['%chk=MBO_' int2str(suffix_L) '_' int2str(suffix_R) '.chk']); %MBO
%         file_contents_current=regexprep(file_contents,'%chk=MBN.chk',['%chk=MBN_' int2str(suffix_L) '_' int2str(suffix_R) '.chk']); %MBN
        file_contents_current=regexprep(file_contents_current,line_old_L,line_new_L);
        file_contents_current=regexprep(file_contents_current,line_old_R,line_new_R);
        
        file_dest_gjf=...
            [path_dest '\' filename_prefix_dest int2str(suffix_L) '_' int2str(suffix_R) file_ext_gjf];
        
%         file_dest_chk=...
%             [path_dest '\' filename_prefix_dest int2str(suffix_L) '_' int2str(suffix_R) file_ext_chk];
        
        fout=fopen(file_dest_gjf,'w+');
        
        fprintf(fout,'%s',file_contents_current);
        
        fclose('all');
        
%         copyfile(file_source_chk,file_dest_chk);
    end
end

fclose('all');