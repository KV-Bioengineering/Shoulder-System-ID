%% dataPrep - function to determine file name, experiment (contraction) type and import corresponding lvm data
function [InAngle,OutTrq,filenm,exptype] = dataPrep(file,path)

    filePath = cat(2,path,file); % construct filepath from file name and dir path

    lvmData = lvm_import(filePath,1); %import lvm data
    
    filenm = extractBefore(file,".");  %Get Name of file before extension
    nmparts = regexp(filenm,'[_()]','split');   %Splits the filename into parts below

    ctnLvl = string(nmparts(1)); %Target Contraction of Trial (%age)
    exptype = string(nmparts(2));  %Trial Contraction Type
    idno = nmparts(3); %Participant ID
    side = string(nmparts(4)); %Left or right arm

    switch exptype
        case {'Ab','Ad'}
            InAngle = lvmData.Segment1.data(:,4);
            OutTrq = lvmData.Segment1.data(:,5);
        case {'In','Ex','IE'}
            InAngle = lvmData.Segment1.data(:,1);
            OutTrq = lvmData.Segment1.data(:,2);
    end
end