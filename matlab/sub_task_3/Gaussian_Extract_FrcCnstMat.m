% Suppressed warnings
%#ok<*SAGROW>

clear all
clc

dim_3_L=0;
dim_3_R=0;
rd_angle_L_MBN=[];
rd_angle_R_MBN=[];
frc_cnst_L_MBN=[];
frc_cnst_R_MBN=[];
molecule_name='MBN';

path_source=['C:\Users\Thomas\Desktop\GaussView 3.09\Myoglobin 2016\' molecule_name '_F_Output\2-29-16'];
path_dest=['C:\Users\Thomas\Desktop\GaussView 3.09\Myoglobin 2016\' molecule_name '_F_Output\2-29-16'];

filename_prefix_source_L_data=[molecule_name '_L_'];
filename_prefix_source_R_data=[molecule_name '_R_'];
file_ext_log='.log';
file_ext_txt='.txt';
file_dest_data=[path_dest '\' filename_prefix_source_L_data file_ext_txt];

line_start_reading=' Force constants in internal coordinates: ';
line_stop_reading=' Final forces over variables, ';

ln_type_1='\s+(?<row>\d+)\s+(?<element_1>[-]*\d+[.]\d+[e]+[+-]*\d+)';
ln_type_2='\s+(?<row>\d+)\s+(?<element_1>[-]*\d+[.]\d+[e]+[+-]*\d+)\s+(?<element_2>[-]*\d+[.]\d+[e]+[+-]*\d+)';
ln_type_3='\s+(?<row>\d+)\s+(?<element_1>[-]*\d+[.]\d+[e]+[+-]*\d+)\s+(?<element_2>[-]*\d+[.]\d+[e]+[+-]*\d+)\s+(?<element_3>[-]*\d+[.]\d+[e]+[+-]*\d+)';
ln_type_4='\s+(?<row>\d+)\s+(?<element_1>[-]*\d+[.]\d+[e]+[+-]*\d+)\s+(?<element_2>[-]*\d+[.]\d+[e]+[+-]*\d+)\s+(?<element_3>[-]*\d+[.]\d+[e]+[+-]*\d+)\s+(?<element_4>[-]*\d+[.]\d+[e]+[+-]*\d+)';
ln_type_5='\s+(?<row>\d+)\s+(?<element_1>[-]*\d+[.]\d+[e]+[+-]*\d+)\s+(?<element_2>[-]*\d+[.]\d+[e]+[+-]*\d+)\s+(?<element_3>[-]*\d+[.]\d+[e]+[+-]*\d+)\s+(?<element_4>[-]*\d+[.]\d+[e]+[+-]*\d+)\s+(?<element_5>[-]*\d+[.]\d+[e]+[+-]*\d+)';

for suffix=-90:5:90
    file_source_L_data=[path_source '\' filename_prefix_source_L_data int2str(suffix) file_ext_log];
    file_source_R_data=[path_source '\' filename_prefix_source_R_data int2str(suffix) file_ext_log];
    
    frc_cnst_mat=[];
    read_flag=0;
    
    f_in=fopen(file_source_L_data);
    f_in_line=fgetl(f_in);
    while ischar(f_in_line)
        if regexp(f_in_line,line_start_reading),read_flag=1;end
        
        if read_flag
            ln_len=length(f_in_line);
            f_in_line=regexprep(f_in_line,'D','e');
            switch ln_len
                case 21;
                    temp_values=regexp(f_in_line,ln_type_1,'names');
                    row=str2num(temp_values.row);
                    col=row;
                    col_range=col:col+4;
                    element_1=str2double(temp_values.element_1);
                    frc_cnst_mat(row,col_range)=[element_1 NaN NaN NaN NaN];
                    element_1=[];
                    temp_values=[];
                case 35;
                    temp_values=regexp(f_in_line,ln_type_2,'names');
                    row=str2num(temp_values.row);
                    element_1=str2double(temp_values.element_1);
                    element_2=str2double(temp_values.element_2);
                    frc_cnst_mat(row,col_range)=[element_1 element_2 NaN NaN NaN];
                    element_1=[];
                    element_2=[];
                    temp_values=[];
                case 49;
                    temp_values=regexp(f_in_line,ln_type_3,'names');
                    row=str2num(temp_values.row);
                    element_1=str2double(temp_values.element_1);
                    element_2=str2double(temp_values.element_2);
                    element_3=str2double(temp_values.element_3);
                    frc_cnst_mat(row,col_range)=[element_1 element_2 element_3 NaN NaN];
                    element_1=[];
                    element_2=[];
                    element_3=[];
                    temp_values=[];
                case 63;
                    temp_values=regexp(f_in_line,ln_type_4,'names');
                    row=str2num(temp_values.row);
                    element_1=str2double(temp_values.element_1);
                    element_2=str2double(temp_values.element_2);
                    element_3=str2double(temp_values.element_3);
                    element_4=str2double(temp_values.element_4);
                    frc_cnst_mat(row,col_range)=[element_1 element_2 element_3 element_4 NaN];
                    element_1=[];
                    element_2=[];
                    element_3=[];
                    element_4=[];
                    temp_values=[];
                case 77;
                    temp_values=regexp(f_in_line,ln_type_5,'names');
                    row=str2num(temp_values.row);
                    element_1=str2double(temp_values.element_1);
                    element_2=str2double(temp_values.element_2);
                    element_3=str2double(temp_values.element_3);
                    element_4=str2double(temp_values.element_4);
                    element_5=str2double(temp_values.element_5);
                    frc_cnst_mat(row,col_range)=[element_1 element_2 element_3 element_4 element_5];
                    element_1=[];
                    element_2=[];
                    element_3=[];
                    element_4=[];
                    element_5=[];
                    temp_values=[];
                otherwise
                    row=[];
                    col=[];
                    col_range=[];
            end
            if regexp(f_in_line, line_stop_reading),read_flag=0;end
        end
        
        f_in_line=fgetl(f_in);
    end
    fclose('all');
    if ~isempty(frc_cnst_mat)
        rd_angle_L_MBN=[rd_angle_L_MBN suffix];
        dim_3_L=dim_3_L+1;
        frc_cnst_L_MBN(:,:,dim_3_L)=frc_cnst_mat;
    end
    
    frc_cnst_mat=[];
    read_flag=0;

    f_in=fopen(file_source_R_data);
    f_in_line=fgetl(f_in);
    while ischar(f_in_line)
        if regexp(f_in_line,line_start_reading),read_flag=1;end
        
        if read_flag
            ln_len=length(f_in_line);
            f_in_line=regexprep(f_in_line,'D','e');
            switch ln_len
                case 21;
                    temp_values=regexp(f_in_line,ln_type_1,'names');
                    row=str2num(temp_values.row);
                    col=row;
                    col_range=col:col+4;
                    element_1=str2double(temp_values.element_1);
                    frc_cnst_mat(row,col_range)=[element_1 NaN NaN NaN NaN];
                    element_1=[];
                    temp_values=[];
                case 35;
                    temp_values=regexp(f_in_line,ln_type_2,'names');
                    row=str2num(temp_values.row);
                    element_1=str2double(temp_values.element_1);
                    element_2=str2double(temp_values.element_2);
                    frc_cnst_mat(row,col_range)=[element_1 element_2 NaN NaN NaN];
                    element_1=[];
                    element_2=[];
                    temp_values=[];
                case 49;
                    temp_values=regexp(f_in_line,ln_type_3,'names');
                    row=str2num(temp_values.row);
                    element_1=str2double(temp_values.element_1);
                    element_2=str2double(temp_values.element_2);
                    element_3=str2double(temp_values.element_3);
                    frc_cnst_mat(row,col_range)=[element_1 element_2 element_3 NaN NaN];
                    element_1=[];
                    element_2=[];
                    element_3=[];
                    temp_values=[];
                case 63;
                    temp_values=regexp(f_in_line,ln_type_4,'names');
                    row=str2num(temp_values.row);
                    element_1=str2double(temp_values.element_1);
                    element_2=str2double(temp_values.element_2);
                    element_3=str2double(temp_values.element_3);
                    element_4=str2double(temp_values.element_4);
                    frc_cnst_mat(row,col_range)=[element_1 element_2 element_3 element_4 NaN];
                    element_1=[];
                    element_2=[];
                    element_3=[];
                    element_4=[];
                    temp_values=[];
                case 77;
                    temp_values=regexp(f_in_line,ln_type_5,'names');
                    row=str2num(temp_values.row);
                    element_1=str2double(temp_values.element_1);
                    element_2=str2double(temp_values.element_2);
                    element_3=str2double(temp_values.element_3);
                    element_4=str2double(temp_values.element_4);
                    element_5=str2double(temp_values.element_5);
                    frc_cnst_mat(row,col_range)=[element_1 element_2 element_3 element_4 element_5];
                    element_1=[];
                    element_2=[];
                    element_3=[];
                    element_4=[];
                    element_5=[];
                    temp_values=[];
                otherwise
                    row=[];
                    col=[];
                    col_range=[];
            end
            if regexp(f_in_line, line_stop_reading),read_flag=0;end
        end
        
        f_in_line=fgetl(f_in);
    end
    fclose('all');
    if ~isempty(frc_cnst_mat)
        rd_angle_R_MBN=[rd_angle_R_MBN suffix];
        dim_3_R=dim_3_R+1;
        frc_cnst_R_MBN(:,:,dim_3_R)=frc_cnst_mat;
    end

end



