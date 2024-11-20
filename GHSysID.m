%% GHSysID - Main function for GH system ID
clear all
% function GHSysID()

    addpath(genpath("Functions"));

    Fs = 1000; % Sample Freq
    Dwin = 3000; %window size for device tf calculation
    Dov = 0.6; % overlap for device tf calculation

    %   Device Data ----------
    [DAA_file,DAA_path] = uigetfile('*.lvm','Select Device AB/AD LVM file');
    % DAA_file = 'Device_Ab(0000).lvm'; % for testing
    % DAA_path = '..\PData\0000\'; % for testing

    [DIE_file,DIE_path] = uigetfile('*.lvm','Select Device Int/Ext LVM file');
    % DIE_file = 'Device_In(0000).lvm'; % for testing
    % DIE_path = '..\PData\0000\'; % for testing

    % 1st column Angle Data, 2nd column Trq Data
    [Device_AA(:,1), Device_AA(:,2)] = dataPrep(DAA_file, DAA_path);
    sFRF_DAA = tfestimate(detrend(Device_AA(:,1),1),detrend(Device_AA(:,2),1),Dwin,Dov*Dwin,[],Fs);

    [Device_IE(:,1), Device_IE(:,2)] = dataPrep(DIE_file, DIE_path);
    sFRF_DIE = tfestimate(detrend(Device_IE(:,1),1),detrend(Device_IE(:,2),1),Dwin,Dov*Dwin,[],Fs);
    
    %%   User Data -----------
    [file,path] = uigetfile('*.lvm','Select LVM file','MultiSelect','on');
    % file = '0%_Ab(0000)L.lvm'; % for testing
    % path = '..\PData\0000\'; % for testing
    
    [InAngle,OutTrq,fileName,exptype] = dataPrep(file,path);

    %% Calculation and Output
    switch exptype
        case {'Ab','Ad'}
            [Pars,mnTrq,VAF_Par,sFRF_par,sFRF,coh,ftfArray] = PCalc(Fs,InAngle,OutTrq,sFRF_DAA);
        case {'In','Ex'}
            [Pars,mnTrq,VAF_Par,sFRF_par,sFRF,coh,ftfArray] = PCalc(Fs,InAngle,OutTrq,sFRF_DIE);
    end

    outpath = string(path)+"output\"+string(fileName)+"\";
    % outpath = "output\"+string(fileName);

    if ~exist(outpath,'dir')
        mkdir(outpath);
    end

    writematrix(Pars,outpath+string(fileName)+"_Par.csv");
    writematrix(mnTrq,outpath+string(fileName)+"_meanTrq.csv");
    writematrix(VAF_Par,outpath+string(fileName)+"_VAFPar.csv");
    writematrix(sFRF_par,outpath+string(fileName)+"sFRFPar.csv");
    writematrix(sFRF,outpath+string(fileName)+"_sFRF.csv");
    writematrix(coh,outpath+string(fileName)+"_coh.csv");
    writematrix(ftfArray,outpath+string(fileName)+"_ftf.csv");

% end