NFiles=5;
wholeN=0;
for ii=1:NFiles
    decN=(ii-1)*5;
    if decN==0
        wholeN=wholeN+1;
    end
    In='';
    Out='NormalDataMBO';
    if decN==0||decN==5
        sampleName=strcat('MBO_',int2str(wholeN),'_0',int2str(decN));
    else
        sampleName=strcat('MBO_',int2str(wholeN),'_',int2str(decN));
    end
    combData(In,Out,sampleName);
end
