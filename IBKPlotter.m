% IBKplotter - Plots IBK vs mean contraction


% path = uigetdir;    %ask use to select folder
% path = "..\Test Data\Yahya Participants\6456"; %for testing
path = "..\PData\0000\output v1.1.2";

%%For making Participant Database (pid_db) structure
pathtype = path+"\*.csv";
% pathtype = path+"\*.txt";   %for testing
files = dir(pathtype);  %gets file list from folder

pararrayLAA = [];
pararrayRAA = [];
pararrayLIE = [];
pararrayRIE = [];
trqarrayLAA = [];
trqarrayRAA = [];
trqarrayLIE = [];
trqarrayRIE = [];

%Parsing through files
for indx=1:length(files)

    filenm = extractBefore(files(indx).name,".");  %Get Name of file before extension
    nmparts = regexp(filenm,'[_()]','split');   %Splits the filename into parts below

    ctnLvl = string(nmparts(1)); %Target Contraction of Trial (%age)
    exptype = string(nmparts(2));  %Trial Contraction Type
    idno = nmparts(3); %Participant ID
    side = string(nmparts(4)); %Left or right arm
    dtype = nmparts(5); % data contained in the file


    % do same thing for mn trq
    switch exptype
        case {'Ab','Ad'}
            % if name contains pars -> logic 1. for all true get data and arrange into pararray
            if ismember(dtype,'Par')
                if ismember(side,'L')
                    pararrayLAA = [pararrayLAA readmatrix(string(files(indx).folder)+"\"+string(files(indx).name))'];
                else
                    pararrayRAA = [pararrayRAA readmatrix(string(files(indx).folder)+"\"+string(files(indx).name))'];
                    % separate L & R
                end
            end

            if ismember(dtype,'meanTrq')
                if ismember(side,'L')
                    trqarrayLAA = [trqarrayLAA readmatrix(string(files(indx).folder)+"\"+string(files(indx).name))];
                else
                    trqarrayRAA = [trqarrayRAA readmatrix(string(files(indx).folder)+"\"+string(files(indx).name))];
                    % separate L & R
                end
            end
        case {'In','Ex'}
            if ismember(dtype,'Par')
                if ismember(side,'L')
                    pararrayLIE = [pararrayLIE readmatrix(string(files(indx).folder)+"\"+string(files(indx).name))'];
                else
                    pararrayRIE = [pararrayRIE readmatrix(string(files(indx).folder)+"\"+string(files(indx).name))'];
                    % separate L & R
                end
            end

            if ismember(dtype,'meanTrq')
                if ismember(side,'L')
                    trqarrayLIE = [trqarrayLIE readmatrix(string(files(indx).folder)+"\"+string(files(indx).name))];
                else
                    trqarrayRIE = [trqarrayRIE readmatrix(string(files(indx).folder)+"\"+string(files(indx).name))];
                    % separate L & R
                end
            end
    end

end

% plot graph

% Inertia vs Trq
IBKfig = figure;
IBKplt = tiledlayout(3,2);

I_IEplt = nexttile([1 1]);
hold on
scatter(trqarrayLIE,pararrayLIE(1,:),[],"blue","DisplayName","Left arm");
scatter(trqarrayRIE,pararrayRIE(1,:),[],"red","DisplayName","Right arm");
title(I_IEplt,"Internal/External Rotation")
ylabel(I_IEplt,"I (kg.m^2)")
hold off

nexttile();
hold on
scatter(trqarrayLAA,pararrayLAA(1,:),[],"blue","DisplayName","Left arm");
scatter(trqarrayRAA,pararrayRAA(1,:),[],"red","DisplayName","Right arm");
title("Abduction/Adduction")
hold off

% Viscous Resistance vs Trq
nexttile
hold on
scatter(trqarrayLIE,pararrayLIE(2,:),[],"blue","DisplayName","Left arm");
scatter(trqarrayRIE,pararrayRIE(2,:),[],"red","DisplayName","Right arm");
ylabel("B (Nm.s/rad)")
hold off

nexttile
hold on
scatter(trqarrayLAA,pararrayLAA(2,:),[],"blue","DisplayName","Left arm");
scatter(trqarrayRAA,pararrayRAA(2,:),[],"red","DisplayName","Right arm");
legend("Location","best")
hold off

% Stiffness vs Trq
nexttile
hold on
scatter(trqarrayLIE,pararrayLIE(3,:),[],"blue","DisplayName","Left arm");
scatter(trqarrayRIE,pararrayRIE(3,:),[],"red","DisplayName","Right arm");

ylabel("K (Nm/rad)")
xlabel("Mean Torque (Nm)")
hold off

nexttile
hold on
scatter(trqarrayLAA,pararrayLAA(3,:),[],"blue","DisplayName","Left arm");
scatter(trqarrayRAA,pararrayRAA(3,:),[],"red","DisplayName","Right arm");
xlabel("Mean Torque (Nm)")
hold off


% % Inertia vs Trq
%         IBKfig = figure;
%         IBKplt = tiledlayout(3,2);
% 
%         I_IEplt = nexttile([1 1]);
%         hold on
%         scatter(trqarrayL,pararrayL(1,:),ParL(1,IntInd),[],"blue","DisplayName","Left arm");
%         scatter(meanTrqL(PIEInd),ParL(1,PIEInd),[],[0.3010 0.7450 0.9330],"filled","DisplayName","Left arm - passive");
% 
%         % scatter([-meanTrqR(ExtInd),-meanTrqR(IntInd)],[ParR(1,ExtInd),ParR(1,IntInd)],[],"red","DisplayName","Right arm - active");
%         % scatter(meanTrqR(PIEInd),ParR(1,PIEInd),[],[0.6350 0.0780 0.1840],"filled","DisplayName","Right arm - passive");
%         title(I_IEplt,"Internal/External Rotation")
%         ylabel(I_IEplt,"I (kg.m^2)")
%         hold off
% 
%         nexttile();
%         hold on
%         scatter([meanTrqL(AdInd),meanTrqL(AbInd)],[ParL(1,AdInd),ParL(1,AbInd)],[],"blue","DisplayName","Left arm - active");
%         scatter(meanTrqL(PAAInd),ParL(1,PAAInd),[],[0.3010 0.7450 0.9330],"filled","DisplayName","Left arm - passive");
% 
%         % scatter([-meanTrqR(AdInd),-meanTrqR(AbInd)],[ParR(1,AbInd),ParR(1,AdInd)],[],"red","DisplayName","Right arm - active");
%         % scatter(meanTrqR(PAAInd),ParR(1,PAAInd),[],[0.6350 0.0780 0.1840],"filled","DisplayName","Right arm - passive");
%         title("Abduction/Adduction")
%         hold off
% 
%         % Viscous Resistance vs Trq
%         nexttile
%         hold on
%         scatter([meanTrqL(ExtInd),meanTrqL(IntInd)],[ParL(3,ExtInd),ParL(3,IntInd)],[],"blue","DisplayName","Left arm - active");
%         scatter(meanTrqL(PIEInd),ParL(3,PIEInd),[],[0.3010 0.7450 0.9330],"filled","DisplayName","Left arm - passive");
% 
%         % scatter([-meanTrqR(ExtInd),-meanTrqR(IntInd)],[ParR(3,ExtInd),ParR(3,IntInd)],[],"red","DisplayName","Right arm - active");
%         % scatter(meanTrqR(PIEInd),ParR(3,PIEInd),[],[0.6350 0.0780 0.1840],"filled","DisplayName","Right arm - passive");
%         % title("Viscous Resistance (Internal/External Rotation)")
%         ylabel("B (Nm.s/rad)")
%         hold off
% 
%         nexttile
%         hold on
%         scatter([meanTrqL(AdInd),meanTrqL(AbInd)],[ParL(2,AdInd),ParL(2,AbInd)],[],"blue","DisplayName","Left arm - active");
%         scatter(meanTrqL(PAAInd),ParL(2,PAAInd),[],[0.3010 0.7450 0.9330],"filled","DisplayName","Left arm - passive");
% 
%         % scatter([-meanTrqR(AdInd),-meanTrqR(AbInd)],[ParR(2,AbInd),ParR(2,AdInd)],[],"red","DisplayName","Right arm - active");
%         % scatter(meanTrqR(PAAInd),ParR(2,PAAInd),[],[0.6350 0.0780 0.1840],"filled","DisplayName","Right arm - passive");
%         legend("Location","best")
%         hold off
% 
%         % Stiffness vs Trq
%         nexttile
%         hold on
%         scatter([meanTrqL(ExtInd),meanTrqL(IntInd)],[ParL(3,ExtInd),ParL(3,IntInd)],[],"blue","DisplayName","Left arm - active");
%         scatter(meanTrqL(PIEInd),ParL(3,PIEInd),[],[0.3010 0.7450 0.9330],"filled","DisplayName","Left arm - passive");
% 
%         % scatter([-meanTrqR(ExtInd),-meanTrqR(IntInd)],[ParR(3,ExtInd),ParR(3,IntInd)],[],"red","DisplayName","Right arm - active");
%         % scatter(meanTrqR(PIEInd),ParR(3,PIEInd),[],[0.6350 0.0780 0.1840],"filled","DisplayName","Right arm - passive");
%         ylabel("K (Nm/rad)")
%         xlabel("Mean Torque (Nm)")
%         hold off
% 
%         nexttile
%         hold on
%         scatter([meanTrqL(AdInd),meanTrqL(AbInd)],[ParL(3,AdInd),ParL(3,AbInd)],[],"blue","DisplayName","Left arm - active");
%         scatter(meanTrqL(PAAInd),ParL(3,PAAInd),[],[0.3010 0.7450 0.9330],"filled","DisplayName","Left arm - passive");
% 
%         % scatter([-meanTrqR(AdInd),-meanTrqR(AbInd)],[ParR(3,AbInd),ParR(3,AdInd)],[],"red","DisplayName","Right arm - active");
%         % scatter(meanTrqR(PAAInd),ParR(3,PAAInd),[],[0.6350 0.0780 0.1840],"filled","DisplayName","Right arm - passive");
%         xlabel("Mean Torque (Nm)")
%         hold off
