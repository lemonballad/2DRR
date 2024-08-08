%#ok<*AGROW> % Suppresses growing array warning
path_source='C:\Users\Thomas\Desktop\GaussView 3.09\Myoglobin 2016\MBO_F_Output\2-29-16';
path_dest='C:\Users\Thomas\Desktop\GaussView 3.09\Myoglobin 2016\MBO_F_Output\2-29-16';
molecule_name='MBO';
filename_prefix_source_L_data=[molecule_name '_L_'];
filename_prefix_source_R_data=[molecule_name '_R_'];
filename_prefix_dest_data=molecule_name;
file_ext_log='.log';
file_ext_txt='.txt';
file_dest_data=[path_dest '\' filename_prefix_dest_data file_ext_txt];

line_start_reading=' Full mass-weighted force constant matrix:';
line_stop_reading=' - Thermochemistry -';
line_1=' Frequencies --\s+(?<w_1>[-]*\d+[.]\d+)\s+(?<w_2>[-]*\d+[.]\d+)\s+(?<w_3>[-]*\d+[.]\d+)';
line_2=' Frc consts  --\s+(?<k_1>[-]*\d+[.]\d+)\s+(?<k_2>[-]*\d+[.]\d+)\s+(?<k_3>[-]*\d+[.]\d+)';
line_3=' Red. masses --\s+(?<rm_1>[-]*\d+[.]\d+)\s+(?<rm_2>[-]*\d+[.]\d+)\s+(?<rm_3>[-]*\d+[.]\d+)';

ws=[];
ks=[];
rm=[];
freq_L_mat=[];
freq_R_mat=[];
force_cnst_L_mat=[];
force_cnst_R_mat=[];
red_mass_L_mat=[];
red_mass_R_mat=[];
rd_angle_L=[];
rd_angle_R=[];

for suffix=-90:5:90 % Loop over files to read
    file_source_L_data=[path_source '\' filename_prefix_source_L_data...
        int2str(suffix) file_ext_log];
    file_source_R_data=[path_source '\' filename_prefix_source_R_data...
        int2str(suffix) file_ext_log];
    
    f_in_L=fopen(file_source_L_data);
    f_in_R=fopen(file_source_R_data);
    
    f_in_L_line=fgetl(f_in_L);
    f_in_R_line=fgetl(f_in_R);
    
    read_flag=0;
    freq=[];
    force_cnst=[];
    red_mass=[];
    while ischar(f_in_L_line) % Loop over lines in file
        if regexp(f_in_L_line,line_start_reading),read_flag=1;end
        
        if read_flag % Begin looking through data
            ws=regexp(f_in_L_line,line_1,'names');
            ks=regexp(f_in_L_line,line_2,'names');
            rm=regexp(f_in_L_line,line_3,'names');
            
            if size(ws,1)>0
                freq=[freq str2double(ws.w_1) str2double(ws.w_2) str2double(ws.w_3)];
                ws=[];
            end
            
            if size(ks,1)>0
                force_cnst=[force_cnst str2double(ks.k_1) str2double(ks.k_2) str2double(ks.k_3)];
                ks=[];
            end
            
            if size(rm,1)>0
                red_mass=[red_mass str2double(rm.rm_1) str2double(rm.rm_2) str2double(rm.rm_3)];
                rm=[];
            end
            
            if regexp(f_in_L_line, line_stop_reading),read_flag=0;end
        end % All data has been collected from this file
        f_in_L_line=fgetl(f_in_L);
    end % loop over lines in file L
    if size(freq)>0,rd_angle_L=[rd_angle_L suffix];end
    freq_L_mat=[freq_L_mat;freq];
    force_cnst_L_mat=[force_cnst_L_mat;force_cnst];
    red_mass_L_mat=[red_mass_L_mat;red_mass];
    
    read_flag=0;
    freq=[];
    force_cnst=[];
    red_mass=[];
    while ischar(f_in_R_line) % Loop over lines in file
        if regexp(f_in_R_line,line_start_reading),read_flag=1;end
        
        if read_flag % Begin looking through data
            ws=regexp(f_in_R_line,line_1,'names');
            ks=regexp(f_in_R_line,line_2,'names');
            rm=regexp(f_in_R_line,line_3,'names');
            
            if size(ws,1)>0
                freq=[freq str2double(ws.w_1) str2double(ws.w_2) str2double(ws.w_3)];
                ws=[];
            end
            
            if size(ks,1)>0
                force_cnst=[force_cnst str2double(ks.k_1) str2double(ks.k_2) str2double(ks.k_3)];
                ks=[];
            end
            
            if size(rm,1)>0
                red_mass=[red_mass str2double(rm.rm_1) str2double(rm.rm_2) str2double(rm.rm_3)];
                rm=[];
            end
            
            if regexp(f_in_R_line, line_stop_reading),read_flag=0;end
        end % All data has been collected from this file
        f_in_R_line=fgetl(f_in_R);
    end % loop over lines in file R
    fclose('all');
    if size(freq)>0,rd_angle_R=[rd_angle_R suffix];end
    freq_R_mat=[freq_R_mat;freq];
    force_cnst_R_mat=[force_cnst_R_mat;force_cnst];
    red_mass_R_mat=[red_mass_R_mat;red_mass];
end % Loop over files to read

% f_out=fopen(file_dest_data,'w+');
%     fprintf(fout_L,'%s',file_contents_current_L);
fclose('all');

rd_angle_L_MBO=rd_angle_L;
freq_L_MBO_mat=freq_L_mat;
force_cnst_L_MBO_mat=force_cnst_L_mat;
red_mass_L_MBO_mat=red_mass_L_mat;

rd_angle_R_MBO=rd_angle_R;
freq_R_MBO_mat=freq_R_mat;
force_cnst_R_MBO_mat=force_cnst_R_mat;
red_mass_R_MBO_mat=red_mass_R_mat;

if 0
freq_MBO_axis=repmat(freq_L_mat(find(rd_angle_L_MBO==0))',size(freq_L_mat,1),1);

force_cnst_L_mat(~isfinite(force_cnst_L_mat))=NaN;
Good_indices=~isnan(force_cnst_L_mat);

FC_L_MBO=griddata(freq_L_mat(Good_indices),...
                            rd_angle_L_MBO(Good_indices),...
                            force_cnst_L_mat(Good_indices),...
                            freq_MBO_axis,rd_angle_L_MBO);
end
