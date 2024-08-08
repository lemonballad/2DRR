function combData(In,Out,sampleName)
fIn=fopen(strcat(In,'.log'));                                                                             % Add sufix to input file
fOut=fopen(strcat(Out,'.mat'));                                                                        % Add sufix to output file
freq=[];
while ~feof(fIn)                                                                                               % Check that input file is not at the end
    readline=fgetl(fIn);                                                                                      % Read line from input.log
    if ~isempty(strfind(readline,' Frequencies --'))
        freq=[freq textscan(readline,...
            strcat(' Frequencies --     ','%f %f %f'))];                                             % Read frequency data
    elseif ~isempty(strfind(readline,' Red. masses --'))
        redM=[redM,textscan(readline,...
            strcat(' Red. masses --     ','%f %f %f'))];                                            % Read reduced mass data
    elseif ~isempty(strfind(readline,' Frc consts  --'))
        frcConsts=[frcConsts,textscan(readline,...
            strcat(' Frc consts  --      ','%f %f %f'))];                                              % Read force constant data
    elseif ~isempty(strfind(readline,' IR Inten    --'))
        IR=[IR,textscan(readline,...
            strcat(' IR Inten    --      ','%f %f %f'))];                                                % Read IR intensity data
    elseif ~isempty(strfind(readline,' Raman Activ --'))
        RamAct=[RamAct,textscan(readline,...
            strcat(' Raman Activ --      ','%f %f %f'))];                                           % Read raman activity data
    elseif ~isempty(strfind(readline,' Depolar (P) --'))
        DePolP=[DePolP,textscan(readline,...
            strcat(' Depolar (P) --      ','%f %f %f'))];                                             % Read depolar data
    elseif ~isempty(strfind(readline,' Depolar (U) --'))        
        DePolU=[DePolU,textscan(readline,...
            strcat(' Depolar (U) --      ','%f %f %f'))];                                             % Read depolar data
    end
end     % End of while
data=genvarname(sampleName);                                                                   % Create variable name
eval([data '=GauData(freq,redM,frcConsts,IR,RamAct,DePolP,DePolU']);    % Assign read in data to structure
save(strcat(Out,'.mat'),'data','-append');                                                       % Save veariables to file
fclose(fIn);                                                                                                       % Close input file
fclose(fOut);                                                                                                    % Close output file
return      % End of Function
