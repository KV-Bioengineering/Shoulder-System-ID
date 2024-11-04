% GHSysID %

% function GHSysID()

    Fs = 1000; % Sample Freq
    % Dwin = 3000; %window size for device tf calculation
    % Dov = 0.8 % overlap for device tf calculation

    %   Device Data ----------
    [DAA_file,DAA_path] = uigetfile('*.lvm','Select Device AB/AD LVM file');

    [DIE_file,DIE_path] = uigetfile('*.lvm','Select Device Int/Ext LVM file');

    % 1st column Angle Data, 2nd column Trq Data
    [Device_AA(:,1), Device_AA(:,2)] = dataPrep(DAA_file, DAA_path);

    [Device_IE(:,1), Device_IE(:,2)] = dataPrep(DIE_file, DIE_path);
    
    %%   User Data -----------
    [file,path] = uigetfile('*.lvm','Select LVM file','MultiSelect','on');
    
    [InAngle,OutTrq, exptype] = dataPrep(file,path);
%%
    switch exptype
        case {'Ab','Ad'}
            % sFRF_D = tfestimate(detrend(Device_AA(:,1),0),detrend(Device_AA(:,2),0),[],[],size(Device_AA,1),Fs);
            [Pararray,sFRF_array,sFRF_parArray,cohArray,VAF_NParray,VAF_ParArray,ftfArray] = ParameterCalc(Fs,InAngle,OutTrq,Device_AA);
        case {'In','Ex'}
            % sFRF_D = tfestimate(detrend(Device_IE(:,1),0),detrend(Device_IE(:,2),0),[],[],size(Device_IE,1),Fs);
            [Pararray,sFRF_array,sFRF_parArray,cohArray,VAF_NParray,VAF_ParArray,ftfArray] = ParameterCalc(Fs,InAngle,OutTrq,Device_IE);
    end
% end