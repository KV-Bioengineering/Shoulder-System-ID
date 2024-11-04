function [InAngle, OutTrq, exptype] = dataPrep(file,path)

    % [file,path] = uigetfile('*.lvm','Select LVM file','MultiSelect','on'); % For testing, this is passed through from main function

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
            InAngle = lvmData.Segment1.data(1:40001,4);
            OutTrq = lvmData.Segment1.data(1:40001,5);
        case {'In','Ex'}
            InAngle = lvmData.Segment1.data(1:40001,1);
            OutTrq = lvmData.Segment1.data(1:40001,2);
    end

    % varargout(1) = exptype;
    
end