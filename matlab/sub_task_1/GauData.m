classdef GauData < handle
    properties
        freq;
        redM;
        frcConsts;
        IR;
        RamAct;
        DePolP;
        DePolU;
    end
    methods
        function data=GauData(freq,redM,frcConsts,IR,RamAct,DePolP,DePolU)
            data.freq=double(freq);
            data.redM=double(redM);
            data.frcConsts=double(frcConsts);
            data.IR=double(IR);
            data.RamAct=double(RamAct);
            data.DePolP=double(DePolP);
            data.DePolU=double(DePolU);
        end
    end
end % End class definition