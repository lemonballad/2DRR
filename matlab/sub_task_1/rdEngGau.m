function [eWN,dis]=rdEngGau(list,pattern1,pattern2)

% 1) Reads a list of file names of Gaussian .log
% 2) Reads energies in hartrees from .log
% 3) Converts energies to cm^-1

first=1;                                                                    % Initialize first iteration flag to true
fail=' Convergence failure -- run terminated.';
index=0;                                                                    % Initialize data index to 0
fileID=fopen(list);                                                         % Open list file
while ~feof(fileID)                                                         % Check if end of file
    tempLine=fgetl(fileID);                                                 % Read in line from list file
    if first==1     % First loop through
        numfiles=sscanf(tempLine,'%i',1);                                   % Read in number of .log files to read
        ehartree(1:numfiles)=double(0);                                     % Initialize energy array to 0s
        dis(1:numfiles)=double(0);
        first=0;                                                            % Set first loop flag to finished
    else            % Subsequent loops through
        tempID=fopen(tempLine);                                             % Open .log file
        index=index+1;                                                      % Increase index for new data
        while ~feof(tempID)                                                 % Check if end of file
            readline=fgetl(tempID);                                         % Read in line from .log file
            if ~isempty(uint8(strfind(readline,pattern1)))                   % Check if line has the desired data
                tempPat1=strcat(pattern1,'%s');
                tempNum1=...
                    str2num(sscanf(readline,tempPat1));                     % Read energy data
                ehartree(index)=tempNum1;
            end
            if ~isempty(uint8(strfind(readline,pattern2)))
                tempPat2=strcat(pattern2,'%f');
                tempNum2=...
                    sscanf(readline,tempPat2);                              % Read distance data
                dis(index)=tempNum2;
            end     % End of if black
            if ~isempty(uint8(strfind(readline,fail)))
                dis(index)=0;
                ehartree(index)=0;
            end
        end     % End Of Gaussian File
        fclose(tempID);
    end     % End if block
end     % End Of List File
fclose(fileID);

eWN=ehartree;%*219474.6313705;                                                % Convert energy in hartrees to wavenumbers
return

