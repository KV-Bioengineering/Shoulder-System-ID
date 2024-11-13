% GHSysID %
clear all
% function GHSysID()

    Fs = 1000; % Sample Freq
    % Dwin = 3000; %window size for device tf calculation
    % Dov = 0.8 % overlap for device tf calculation

    %   Device Data ----------
    [DAA_file,DAA_path] = uigetfile('*.lvm','Select Device AB/AD LVM file');
    % DAA_file = 'Device_Ab(0000).lvm'; % for testing
    % DAA_path = '..\PData\0000\'; % for testing

    [DIE_file,DIE_path] = uigetfile('*.lvm','Select Device Int/Ext LVM file');
    % DIE_file = 'Device_In(0000).lvm'; % for testing
    % DIE_path = '..\PData\0000\'; % for testing

    % 1st column Angle Data, 2nd column Trq Data
    [Device_AA(:,1), Device_AA(:,2)] = dataPrep(DAA_file, DAA_path);

    [Device_IE(:,1), Device_IE(:,2)] = dataPrep(DIE_file, DIE_path);
    
    %%   User Data -----------
    % [file,path] = uigetfile('*.lvm','Select LVM file','MultiSelect','on');
    file = '0%_Ab(0000)L.lvm';
    path = '..\PData\0000\';
    
    [InAngle,OutTrq,fileName,exptype] = dataPrep(file,path);

    %% Calculation and Output
    switch exptype
        case {'Ab','Ad'}
            % sFRF_D = tfestimate(detrend(Device_AA(:,1),0),detrend(Device_AA(:,2),0),[],[],size(Device_AA,1),Fs);
            % [Pararray,sFRF_array,sFRF_parArray,cohArray,VAF_NParray,VAF_ParArray,ftfArray] = ParameterCalc(Fs,InAngle,OutTrq,Device_AA);
            [Pars,VAF_Par,sFRF_par,sFRF,coh,ftfArray] = PCalc(Fs,InAngle,OutTrq,Device_AA);
        case {'In','Ex'}
            % sFRF_D = tfestimate(detrend(Device_IE(:,1),0),detrend(Device_IE(:,2),0),[],[],size(Device_IE,1),Fs);
            % [Pararray,sFRF_array,sFRF_parArray,cohArray,VAF_NParray,VAF_ParArray,ftfArray] = ParameterCalc(Fs,InAngle,OutTrq,Device_IE);
            [Pars,VAF_Par,sFRF_par,sFRF,coh,ftfArray] = PCalc(Fs,InAngle,OutTrq,Device_IE);
    end

    outpath = string(path)+"output\"+string(fileName);
    % outpath = "output\"+string(fileName);

    writematrix(VAF_Par,outpath+"_VAFPar.xls");
    writematrix(sFRF_par,outpath+"sFRFPar.xls");
    writematrix(sFRF,outpath+"_sFRF.xls");
    writematrix(coh,outpath+"_coh.xls");
    writematrix(ftfArray,outpath+"_ftf.xls");

% end