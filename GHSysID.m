%% GHSysID - Main function for GH system ID
function GHSysID()

    addpath(genpath("Functions"));

    Fs = 1000; % Sample Freq
    Dwin = 3000; %window size for device tf calculation
    Dov = 0.6; % overlap for device tf calculation

    %   Select Output folder
    outpath = uigetdir('','Select Output Directory');
    % outpath = "..\PData\0000\output"; % for testing

    %%  Device Data ----------
    % [DAA_file,DAA_path] = uigetfile('*.lvm','Select Device AB/AD LVM file');
    DAA_file = 'Device_Ab(0000).lvm'; % for testing
    DAA_path = '..\PData\0000\'; % for testing

    % [DIE_file,DIE_path] = uigetfile('*.lvm','Select Device Int/Ext LVM file');
    DIE_file = 'Device_In(0000).lvm'; % for testing
    DIE_path = '..\PData\0000\'; % for testing

    % 1st column Angle Data, 2nd column Trq Data
    [Device_AA(:,1), Device_AA(:,2)] = dataPrep(DAA_file, DAA_path);
    sFRF_DAA = tfestimate(detrend(Device_AA(:,1),1),detrend(Device_AA(:,2),1),Dwin,Dov*Dwin,10001,Fs);

    [Device_IE(:,1), Device_IE(:,2)] = dataPrep(DIE_file, DIE_path);
    sFRF_DIE = tfestimate(detrend(Device_IE(:,1),1),detrend(Device_IE(:,2),1),Dwin,Dov*Dwin,10001,Fs);
    
    %%  User Data -----------
    [file,path] = uigetfile('*.lvm','Select LVM file','MultiSelect','on');
    % file = '0%_Ab(0000)L.lvm'; % for testing
    % path = '..\PData\0000\'; % for testing

    %% Calculation and Output
    % In case a single file is selected
    if ischar(file)
        file = cellstr(file);
    end

    % For loop is for multiple file selections
    for fno = 1:length(file)
        [InAngle,OutTrq,fileName,exptype] = dataPrep(char(file(fno)),path);

        switch exptype
            case {'Ab','Ad'}
                [Pars,mnTrq,VAF_NP,VAF_Par,sFRF_par,sFRF,coh,ftfArray] = PCalc(Fs,InAngle,OutTrq,sFRF_DAA);
            case {'In','Ex'}
                [Pars,mnTrq,VAF_NP,VAF_Par,sFRF_par,sFRF,coh,ftfArray] = PCalc(Fs,InAngle,OutTrq,sFRF_DIE);
        end

        writematrix(Pars,outpath+"\"+string(fileName)+"_Par.csv");
        writematrix(mnTrq,outpath+"\"+string(fileName)+"_meanTrq.csv");
        writematrix(VAF_NP,outpath+"\"+string(fileName)+"_VAFNP.csv");
        writematrix(VAF_Par,outpath+"\"+string(fileName)+"_VAFPar.csv");
        writematrix(sFRF_par,outpath+"\"+string(fileName)+"_sFRFPar.csv");
        writematrix(sFRF,outpath+"\"+string(fileName)+"_sFRF.csv");
        writematrix(coh,outpath+"\"+string(fileName)+"_coh.csv");
        writematrix(ftfArray,outpath+"\"+string(fileName)+"_ftf.csv");
    end

end