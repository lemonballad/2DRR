%%
% Suppress warnings
%#ok<*SAGROW> 
%#ok<*AGROW>

% %% Define paths, FOR MBN
% path_source='C:\Users\Thomas\Desktop\GaussView 3.09\Myoglobin 2016\MBN_Grid_Output\4-6-16';
% molecule='MBN';
%% Define paths, FOR MBO
path_source='C:\Users\Thomas\Desktop\GaussView 3.09\Myoglobin 2016\MBO_Grid_Output\4-6-16';
molecule='MBO';
file_ext_log='.log';

%% Define string expressions to read and flags for reading
% String for reading normal mode cartesian coordinates lines
NM_carts_str='\s+(?<atom>\d+)\s+\d+';
for icart=1:3 % Loop over normal mode per line indices
    for xyz=['x' 'y' 'z'] % Loop over normal mode cartesian coordinate directions
        NM_cart_str_new=['\s+(?<' xyz int2str(icart) '>[-]*\d+[.]+\d+)'];
        NM_carts_str=[NM_carts_str NM_cart_str_new]; % Build full string
    end % End loop over nm cart coords directions
end % End loop over nm/line indices

% String for reading force constant line
frc_cnsts_str=' Frc consts  --\s+(?<k1>[-]*\d+[.]\d+)\s+(?<k2>[-]*\d+[.]\d+)\s+(?<k3>[-]*\d+[.]\d+)';
% String for reading vibrational frequency line
freqs_str=' Frequencies --\s+(?<w1>[-]*\d+[.]\d+)\s+(?<w2>[-]*\d+[.]\d+)\s+(?<w3>[-]*\d+[.]\d+)';
freqs_start_str=' Full mass-weighted force constant matrix:'; % Start saving data here
freqs_end_str=' - Thermochemistry -'; % Stop reading data here
% String for reading mode numbers line
modes_str='\s+(?<m1>\d+)\s+(?<m2>\d+)\s+(?<m3>\d+)';
% String for reading reduced masses line
red_masses_str=' Red. masses --\s+(?<rm1>[-]*\d+[.]\d+)\s+(?<rm2>[-]*\d+[.]\d+)\s+(?<rm3>[-]*\d+[.]\d+)';

%% Read Equilibrium file
freq_read_flag=false; % Should the line be read as possible normal mode data
counter=0;
last_counter=-1;

% Define file to open
file_name=[path_source '\' molecule '_0_0' file_ext_log];
f_in=fopen(file_name); % Open grid file
line=fgetl(f_in); % Read in first line of current grid file
% Loop over lines in file
while ischar(line) % Checking for empty file or end of file
    % Set flag for reading if in NM section of log file
    if regexp(line,freqs_start_str),freq_read_flag=true;end;
    
    if freq_read_flag % If in NM section of log file start collecting data
        mode_data=regexp(line,modes_str,'names'); % Read modes
        % Write normal mode data to matrices
        if length(line)==69 % Does this line conain mode data
            if size(mode_data,1)~=0 % Check to see if mode data was read in
                if counter>last_counter % If there is more recent mode data restart mode count
                    last_counter=counter;
                    modes=[];
                end
                modes=[modes str2double(mode_data.m1) str2double(mode_data.m2) str2double(mode_data.m3)];
            end
        end
        freq_data=regexp(line,freqs_str,'names'); % Read frequencies
        red_mass_data=regexp(line,red_masses_str,'names'); % Read reduced masses
        frc_data=regexp(line,frc_cnsts_str,'names'); % Read force constants
        NM_cart_data=regexp(line,NM_carts_str,'names'); % Read normal mode cartesian coordinates
        if size(NM_cart_data,1)~=0 % Check to see if cartesian data was read
            iatom=str2double(NM_cart_data.atom);
            NM_carts(modes(end)-2,iatom,1:3)=[str2double(NM_cart_data.x1) str2double(NM_cart_data.y1) str2double(NM_cart_data.z1)];
            NM_carts(modes(end)-1,iatom,1:3)=[str2double(NM_cart_data.x2) str2double(NM_cart_data.y2) str2double(NM_cart_data.z2)];
            NM_carts(modes(end),iatom,1:3)=[str2double(NM_cart_data.x3) str2double(NM_cart_data.y3) str2double(NM_cart_data.z3)];
        end
        
        if size(freq_data,1)~=0 % Check to see if frequency data was read
            %                     freq(modes(end)-2:modes(end),iL,iR)=[str2double(freq_data.w1) str2double(freq_data.w2) str2double(freq_data.w3)];
            freq(modes(end-2:end))=[str2double(freq_data.w1) str2double(freq_data.w2) str2double(freq_data.w3)];
        end
        if size(red_mass_data,1)~=0 % Check to see if reduced mass data was read
            %                     red_mass(modes(end)-2:modes(end),iL,iR)=[str2double(red_mass_data.rm1) str2double(red_mass_data.rm2) str2double(red_mass_data.rm3)];
            red_mass(modes(end-2:end))=[str2double(red_mass_data.rm1) str2double(red_mass_data.rm2) str2double(red_mass_data.rm3)];
        end
        if size(frc_data,1)~=0 % Check to see if force constant data was read
            %                     frc_cnst(modes(end)-2:modes(end),iL,iR)=[str2double(frc_data.k1) str2double(frc_data.k2) str2double(frc_data.k3)];
            frc_cnst(modes(end-2:end))=[str2double(frc_data.k1) str2double(frc_data.k2) str2double(frc_data.k3)];
        end
        
        % Reinitialize data structures for reading strings
        mode_data=[];
        freq_data=[];
        red_mass_data=[];
        frc_data=[];
        NM_cart_data=[];
        if regexp(line,freqs_end_str) % Acknowledge all current NM data has been read
            freq_read_flag=false;
            counter=counter+1; % Number of times NM data has been reported in LOG
        end
    end
    line=fgetl(f_in); % Read in next line from file
    
end % End loop over lines of current files
NMC0(:,:)=NM_carts(46,:,:);
normNMC0=sqrt(sum(sum(NMC0.^2)));
NMC0=NMC0/normNMC0;
freq0=freq(46);
frc_cnst0=frc_cnst(46);
nmodes=length(modes);
red_m0=red_mass(46);
clear NM_carts freq frc_cnst red_mass modes;

fclose('all'); % Close all open files. NECESSARY to keep program from eating up memory

%% Read grid files
% Initialize variable
L_prop=-180:5:175;
R_prop=-180:5:175;
DATA_46=zeros(3,length(L_prop),length(R_prop));
frequencies=zeros(nmodes,length(L_prop),length(R_prop));
force_constants=zeros(nmodes,length(L_prop),length(R_prop));
reduced_masses=zeros(nmodes,length(L_prop),length(R_prop));
coordinates=zeros(nmodes,77,3,length(L_prop),length(R_prop));
for iL=1:length(L_prop) % Loop over left propionate dihedral angle
    for iR=1:length(R_prop) % Loop over right propionate dihedral angle
        freq_read_flag=false; % Should the line be read as possible normal mode data
        freq=[];
        frc_cnst=[];
        NM_carts=[];
        red_mass=[];

        % Define file to open
        file_name=[path_source '\' molecule '_' int2str(L_prop(iL)) '_' int2str(R_prop(iR)) file_ext_log];
        f_in=fopen(file_name); % Open grid file
        line=fgetl(f_in); % Read in first line of current grid file
        tic
        % Loop over lines in file
        while ischar(line) % Checking for empty file or end of file
            % Set flag for reading if in NM section of log file
            if regexp(line,freqs_start_str),freq_read_flag=true;end;
            
            if freq_read_flag % If in NM section of log file start collecting data
                mode_data=regexp(line,modes_str,'names'); % Read modes              
                % Write normal mode data to matrices
                if length(line)==69 && size(mode_data,1)~=0 % Does this line contain mode data
                    last_mode=str2double(mode_data.m3);
                end
                
                freq_data=regexp(line,freqs_str,'names'); % Read frequencies
                red_mass_data=regexp(line,red_masses_str,'names'); % Read reduced masses
                frc_data=regexp(line,frc_cnsts_str,'names'); % Read force constants
                NM_cart_data=regexp(line,NM_carts_str,'names'); % Read normal mode cartesian coordinates
                if size(NM_cart_data,1)~=0 % Check to see if cartesian data was read
                    iatom=str2double(NM_cart_data.atom);
                    NM_carts(last_mode-2,iatom,1:3)=[str2double(NM_cart_data.x1) str2double(NM_cart_data.y1) str2double(NM_cart_data.z1)];
                    NM_carts(last_mode-1,iatom,1:3)=[str2double(NM_cart_data.x2) str2double(NM_cart_data.y2) str2double(NM_cart_data.z2)];
                    NM_carts(last_mode,iatom,1:3)=[str2double(NM_cart_data.x3) str2double(NM_cart_data.y3) str2double(NM_cart_data.z3)];
                end
                
                if size(freq_data,1)~=0 % Check to see if frequency data was read
                    freq(last_mode-2:last_mode)=[str2double(freq_data.w1) str2double(freq_data.w2) str2double(freq_data.w3)];
                end
                if size(red_mass_data,1)~=0 % Check to see if reduced mass data was read
                    red_mass(last_mode-2:last_mode)=[str2double(red_mass_data.rm1) str2double(red_mass_data.rm2) str2double(red_mass_data.rm3)];
                end
                if size(frc_data,1)~=0 % Check to see if force constant data was read
                    frc_cnst(last_mode-2:last_mode)=[str2double(frc_data.k1) str2double(frc_data.k2) str2double(frc_data.k3)];
                end
                
                % Reinitialize data structures for reading strings
                mode_data=[];
                freq_data=[];
                red_mass_data=[];
                frc_data=[];
                NM_cart_data=[];
                if regexp(line,freqs_end_str) % Acknowledge all current NM data has been read
                    freq_read_flag=false;
                    counter=counter+1; % Number of times NM data has been reported in LOG
                end
            end
            line=fgetl(f_in); % Read in next line from file
            
        end % End loop over lines of current files
        
        fclose('all'); % Close all open files. NECESSARY to keep program from eating up memory
        
        if ~isempty(freq)
            frequencies(:,iL,iR)=freq;
            force_constants(:,iL,iR)=frc_cnst;
            reduced_masses(:,iL,iR)=red_mass;
            coordinates(1:nmodes,1:77,1:3,iL,iR)=NM_carts;
            max_dot_p=0;
            normNMC=sqrt(sum(sum(NM_carts.^2,3),2));
            normNMC_mat=repmat(normNMC,1,size(NM_carts,2),3);
            NMC=NM_carts./normNMC_mat;
            for imode=35:56
                current_dot_p=sum(sum(NMC0.*permute(NMC(imode,:,:),[2 3 1])));
                if max_dot_p<current_dot_p
                    max_dot_p=current_dot_p;
                    DATA_46(1:5,iL,iR)=[freq(imode) frc_cnst(imode) red_mass(imode) imode current_dot_p];
                end
            end
        else
            DATA_46(1:5,iL,iR)=[NaN NaN NaN NaN NaN];
            frequencies(:,iL,iR)=NaN;
            force_constants(:,iL,iR)=NaN;
            reduced_masses(:,iL,iR)=NaN;
            coordinates(1:nmodes,1:77,1:3,iL,iR)=NaN;
        end
        toc
    end % End loop over right propionate dihedral angle
end % End loop over left propionate dihedral angle



