%% PCalc - function that estimates the parameters of the GH joint, as well as calculates the parametric Variance accounted for
function [Pararray,MnTrq,VAF_NPArray,VAF_ParArray,sFRF_parArray,sFRF_array,cohArray,ftfArray] = PCalc(Fs,inAnglData,outTrqData,sFRF_D)
    Ts = 1/Fs;
    t=(0:10000)*Ts; %time vector
    t=t';

    NSplt = 10001; % size of segments

    %Butterworth Filter initialization
    fc=10; %Cuttoff Freq
    [b,a] = butter(2,fc/(Fs/2),'low');

    %   Angle & Torque Data -----------------
    % Split Data into 4 data sets
    % Angle Data
    inAngle(:,1) = inAnglData(1:10001);
    inAngle(:,2) = inAnglData(10001:20001);
    inAngle(:,3) = inAnglData(20001:30001);
    inAngle(:,4) = inAnglData(30001:40001);

    % Detrended Angle Data
    dtAngle(:,1) = detrend(inAnglData(1:10001),1);
    dtAngle(:,2) = detrend(inAnglData(10001:20001),1);
    dtAngle(:,3) = detrend(inAnglData(20001:30001),1);
    dtAngle(:,4) = detrend(inAnglData(30001:40001),1);

    % Device Torque
    dTrq(:,1) = TrqCalc(sFRF_D,inAngle(:,1));
    dTrq(:,2) = TrqCalc(sFRF_D,inAngle(:,2));
    dTrq(:,3) = TrqCalc(sFRF_D,inAngle(:,3));
    dTrq(:,4) = TrqCalc(sFRF_D,inAngle(:,4));

    %Alternate device torque calculation
    % dTrq(:,1) = TrqConv(sFRF_D,inAngle(:,1));
    % dTrq(:,2) = TrqConv(sFRF_D,inAngle(:,2));
    % dTrq(:,3) = TrqConv(sFRF_D,inAngle(:,3));
    % dTrq(:,4) = TrqConv(sFRF_D,inAngle(:,4));

    % Output Torque (Contraction Torque)
    outTrq(:,1) = outTrqData(1:10001) - dTrq(:,1);
    outTrq(:,2) = outTrqData(10001:20001) - dTrq(:,2);
    outTrq(:,3) = outTrqData(20001:30001) - dTrq(:,3);
    outTrq(:,4) = outTrqData(30001:40001) - dTrq(:,4);

    % Detrended Contraction Torque
    dtTrq(:,1) = detrend(outTrq(:,1),1);
    dtTrq(:,2) = detrend(outTrq(:,2),1);
    dtTrq(:,3) = detrend(outTrq(:,3),1);
    dtTrq(:,4) = detrend(outTrq(:,4),1);

    %%  Transfer Function Estimation -----------------
    win = []; % window size
    ov = [];
    % ov = 0.8*win; % overlap

    % [sFRF,ftf] = tfestimate(inAngl,outTrq,win,ov.*win,f,Fs);

    [sFRF1,ftf1] = tfestimate(dtAngle(:,1),dtTrq(:,1),win,ov,NSplt,Fs); % cross power spectral density between the input and the output
    [sFRF2,ftf2] = tfestimate(dtAngle(:,2),dtTrq(:,2),win,ov,NSplt,Fs);
    [sFRF3,ftf3] = tfestimate(dtAngle(:,3),dtTrq(:,3),win,ov,NSplt,Fs);
    [sFRF4,ftf4] = tfestimate(dtAngle(:,4),dtTrq(:,4),win,ov,NSplt,Fs);

    %   Coherence ------------------------
    coh1 = mscohere(dtAngle(:,1),dtTrq(:,1),win,ov,NSplt,Fs);
    coh2 = mscohere(dtAngle(:,2),dtTrq(:,2),win,ov,NSplt,Fs);
    coh3 = mscohere(dtAngle(:,3),dtTrq(:,3),win,ov,NSplt,Fs);
    coh4 = mscohere(dtAngle(:,4),dtTrq(:,4),win,ov,NSplt,Fs);

    %% Non-Parametric VAF

    fftAngle2s_1 = fft(tukeywin(NSplt,0.02).*dtAngle(:,1)); % 2-sided Angle in freq domain
    fftAngle2s_2 = fft(tukeywin(NSplt,0.02).*dtAngle(:,2)); % 2-sided Angle in freq domain
    fftAngle2s_3 = fft(tukeywin(NSplt,0.02).*dtAngle(:,3)); % 2-sided Angle in freq domain
    fftAngle2s_4 = fft(tukeywin(NSplt,0.02).*dtAngle(:,4)); % 2-sided Angle in freq domain

    % fftAngle = fftAngle2s(1:size(fftAngle2s,1)/2+1);

    fftAngle_1 = fftAngle2s_1(1:size(fftAngle2s_1,1)/2+1); % 1-sided Angle
    fftAngle_2 = fftAngle2s_2(1:size(fftAngle2s_2,1)/2+1); % 1-sided Angle
    fftAngle_3 = fftAngle2s_3(1:size(fftAngle2s_3,1)/2+1); % 1-sided Angle
    fftAngle_4 = fftAngle2s_4(1:size(fftAngle2s_4,1)/2+1); % 1-sided Angle

    % fftTrq = fftAngle.*sFRF';
    
    fftTrq_1 = fftAngle_1.*sFRF1; % calculate Trq in freq domain
    fftTrq_2 = fftAngle_2.*sFRF2;
    fftTrq_3 = fftAngle_3.*sFRF3;
    fftTrq_4 = fftAngle_4.*sFRF4;

    % fftTrq2s = [fftTrq(1); fftTrq(2:end); flipud(conj(fftTrq(2:end)))];

    fftTrq2s_1 = [fftTrq_1(1); fftTrq_1(2:end); flipud(conj(fftTrq_1(2:end)))]; % convert to 2 sided Trq
    fftTrq2s_2 = [fftTrq_2(1); fftTrq_2(2:end); flipud(conj(fftTrq_2(2:end)))];
    fftTrq2s_3 = [fftTrq_3(1); fftTrq_3(2:end); flipud(conj(fftTrq_3(2:end)))];
    fftTrq2s_4 = [fftTrq_4(1); fftTrq_4(2:end); flipud(conj(fftTrq_4(2:end)))];

    % Trq_NP = ifft(fftTrq2s);

    %Take real value of Trq, there seems to be some minor error (order of 10^-15) from sFRF
    % Trq_NP1 = real(ifft(fftTrq2s_1)); % ifft of freq-domain Trq to get time-domain Trq
    % Trq_NP2 = real(ifft(fftTrq2s_2));
    % Trq_NP3 = real(ifft(fftTrq2s_3));
    % Trq_NP4 = real(ifft(fftTrq2s_4));

    Trq_NP1 = ifft(fftTrq2s_1); % ifft of freq-domain Trq to get time-domain Trq
    Trq_NP2 = ifft(fftTrq2s_2);
    Trq_NP3 = ifft(fftTrq2s_3);
    Trq_NP4 = ifft(fftTrq2s_4);

    Trq_filt1 = filtfilt(b,a,dtTrq(:,1));
    Trq_filt2 = filtfilt(b,a,dtTrq(:,2));
    Trq_filt3 = filtfilt(b,a,dtTrq(:,3));
    Trq_filt4 = filtfilt(b,a,dtTrq(:,4));

    Trq_NPFilt1 = filtfilt(b,a,Trq_NP1);
    Trq_NPFilt2 = filtfilt(b,a,Trq_NP2);
    Trq_NPFilt3 = filtfilt(b,a,Trq_NP3);
    Trq_NPFilt4 = filtfilt(b,a,Trq_NP4);

    VAF_NP1 = 100*(1-((var(Trq_filt1-Trq_NPFilt1)))./var(Trq_filt1));
    VAF_NP2 = 100*(1-((var(Trq_filt2-Trq_NPFilt2)))./var(Trq_filt2));
    VAF_NP3 = 100*(1-((var(Trq_filt3-Trq_NPFilt3)))./var(Trq_filt3));
    VAF_NP4 = 100*(1-((var(Trq_filt4-Trq_NPFilt4)))./var(Trq_filt4));

    %%   Parameter Estimation ---------------

    E0=[0.2,2,20];
    opts = optimoptions("lsqnonlin","Algorithm","levenberg-marquardt","StepTolerance",1.000000e-20,"Display","iter","MaxFunctionEvaluations",6000,"MaxIterations",1000);

    %Y-data for parameter estimation
    ftfLim = find(ftf1<10,1,'last');

    w_Splt = 2*pi*ftf1(1:ftfLim);

    ydata1 = sFRF1(1:ftfLim);
    ydata2 = sFRF2(1:ftfLim);
    ydata3 = sFRF3(1:ftfLim);
    ydata4 = sFRF4(1:ftfLim);

    cohydata1 = coh1(1:ftfLim);
    cohydata2 = coh2(1:ftfLim);
    cohydata3 = coh3(1:ftfLim);
    cohydata4 = coh4(1:ftfLim);

    % err1 = @(E1) cohydata1.*(abs(ydata1)+angle(ydata1) - (abs(E1(1)*((i*w_Splt).^2)+E1(2)*i*w_Splt+E1(3))+angle(E1(1)*((i*w_Splt).^2)+E1(2)*i*w_Splt+E1(3))));
    % err2 = @(E2) cohydata2.*(abs(ydata2)+angle(ydata2) - (abs(E2(1)*((i*w_Splt).^2)+E2(2)*i*w_Splt+E2(3))+angle(E2(1)*((i*w_Splt).^2)+E2(2)*i*w_Splt+E2(3))));
    % err3 = @(E3) cohydata3.*(abs(ydata3)+angle(ydata3) - (abs(E3(1)*((i*w_Splt).^2)+E3(2)*i*w_Splt+E3(3))+angle(E3(1)*((i*w_Splt).^2)+E3(2)*i*w_Splt+E3(3))));
    % err4 = @(E4) cohydata4.*(abs(ydata4)+angle(ydata4) - (abs(E4(1)*((i*w_Splt).^2)+E4(2)*i*w_Splt+E4(3))+angle(E4(1)*((i*w_Splt).^2)+E4(2)*i*w_Splt+E4(3))));

    err1 = @(E1) [cohydata1.*(abs(ydata1)-abs((E1(1)*((i*w_Splt).^2)+E1(2)*i*w_Splt+E1(3)))); cohydata1.*(angle(ydata1)-angle((E1(1)*((i*w_Splt).^2)+E1(2)*i*w_Splt+E1(3))))];
    err2 = @(E2) [cohydata2.*(abs(ydata2)-abs((E2(1)*((i*w_Splt).^2)+E2(2)*i*w_Splt+E2(3)))); cohydata2.*(angle(ydata2)-angle((E2(1)*((i*w_Splt).^2)+E2(2)*i*w_Splt+E2(3))))];
    err3 = @(E3) [cohydata3.*(abs(ydata3)-abs((E3(1)*((i*w_Splt).^2)+E3(2)*i*w_Splt+E3(3)))); cohydata3.*(angle(ydata3)-angle((E3(1)*((i*w_Splt).^2)+E3(2)*i*w_Splt+E3(3))))];
    err4 = @(E4) [cohydata4.*(abs(ydata4)-abs((E4(1)*((i*w_Splt).^2)+E4(2)*i*w_Splt+E4(3)))); cohydata4.*(angle(ydata4)-angle((E4(1)*((i*w_Splt).^2)+E4(2)*i*w_Splt+E4(3))))];

    %Least squares fit to estimate IBK
    % E = lsqnonlin(err,E0,[],[],options);
    [E1,resnorm1,residual1] = lsqnonlin(err1,E0,[],[],opts);
    [E2,resnorm2,residual2] = lsqnonlin(err2,E0,[],[],opts);
    [E3,resnorm3,residual3] = lsqnonlin(err3,E0,[],[],opts);
    [E4,resnorm4,residual4] = lsqnonlin(err4,E0,[],[],opts);

    %% Parametric Model and VAF
    fPar = 0:Fs/10001:500; %frequency vector for parametric sFRF, used for graphing

    cIRF1_par = ((2*(-i)*sin(i*(t*(E1(2)^2 - 4*E1(1)*E1(3))^(1/2))/(2*E1(1))).*exp(-(E1(2)*t)/(2*E1(1))))/(E1(2)^2 - 4*E1(1)*E1(3))^(1/2));  %Simulated Compliance IRF (Parametric) <- from the inverse Laplace of the compliance formula
    cIRF2_par = ((2*(-i)*sin(i*(t*(E2(2)^2 - 4*E2(1)*E2(3))^(1/2))/(2*E2(1))).*exp(-(E2(2)*t)/(2*E2(1))))/(E2(2)^2 - 4*E2(1)*E2(3))^(1/2));
    cIRF3_par = ((2*(-i)*sin(i*(t*(E3(2)^2 - 4*E3(1)*E3(3))^(1/2))/(2*E3(1))).*exp(-(E3(2)*t)/(2*E3(1))))/(E3(2)^2 - 4*E3(1)*E3(3))^(1/2));
    cIRF4_par = ((2*(-i)*sin(i*(t*(E4(2)^2 - 4*E4(1)*E4(3))^(1/2))/(2*E4(1))).*exp(-(E4(2)*t)/(2*E4(1))))/(E4(2)^2 - 4*E4(1)*E4(3))^(1/2));

    cFRF1_par=fft(cIRF1_par)/Fs;   % Parametric Compliance FRF
    sFRF1_par=1./cFRF1_par;        % Parametric Stiffness FRF
    sIRF1_parNrot=ifft(sFRF1_par);
    sIRF1_par=ifftshift(sIRF1_parNrot);

    cFRF2_par=fft(cIRF2_par)/Fs;
    sFRF2_par=1./cFRF2_par;
    sIRF2_parNrot=ifft(sFRF2_par);
    sIRF2_par=ifftshift(sIRF2_parNrot);

    cFRF3_par=fft(cIRF3_par)/Fs;
    sFRF3_par=1./cFRF3_par;
    sIRF3_parNrot=ifft(sFRF3_par);
    sIRF3_par=ifftshift(sIRF3_parNrot);

    cFRF4_par=fft(cIRF4_par)/Fs;
    sFRF4_par=1./cFRF4_par;
    sIRF4_parNrot=ifft(sFRF4_par);
    sIRF4_par=ifftshift(sIRF4_parNrot);

    Trq_sim1=conv(tukeywin(NSplt,0.02).*dtAngle(:,1),sIRF1_par,'same'); % Simulated Intrinsic torque
    Trq_sim2=conv(tukeywin(NSplt,0.02).*dtAngle(:,2),sIRF2_par,'same');
    Trq_sim3=conv(tukeywin(NSplt,0.02).*dtAngle(:,3),sIRF3_par,'same');
    Trq_sim4=conv(tukeywin(NSplt,0.02).*dtAngle(:,4),sIRF4_par,'same');

    %   Filtered Torques -----------------

    Trq_simFilt1 = filtfilt(b,a,Trq_sim1);
    Trq_simFilt2 = filtfilt(b,a,Trq_sim2);
    Trq_simFilt3 = filtfilt(b,a,Trq_sim3);
    Trq_simFilt4 = filtfilt(b,a,Trq_sim4);

    %Parametric model Variance Accounted For

    VAF_par1 = 100*(1-((var(Trq_filt1-Trq_simFilt1)))./var(Trq_filt1));
    VAF_par2 = 100*(1-((var(Trq_filt2-Trq_simFilt2)))./var(Trq_filt2));
    VAF_par3 = 100*(1-((var(Trq_filt3-Trq_simFilt3)))./var(Trq_filt3));
    VAF_par4 = 100*(1-((var(Trq_filt4-Trq_simFilt4)))./var(Trq_filt4));

    %% Function output
    % Pararray contains IBK values of each set
    Pararray(1,:) = E1;
    Pararray(2,:) = E2;
    Pararray(3,:) = E3;
    Pararray(4,:) = E4;

    MnTrq(1) = mean(outTrq(:,1));
    MnTrq(2) = mean(outTrq(:,2));
    MnTrq(3) = mean(outTrq(:,3));
    MnTrq(4) = mean(outTrq(:,4));

    VAF_NPArray(1)=VAF_NP1;
    VAF_NPArray(2)=VAF_NP2;
    VAF_NPArray(3)=VAF_NP3;
    VAF_NPArray(4)=VAF_NP4;

    VAF_ParArray(1)=VAF_par1;
    VAF_ParArray(2)=VAF_par2;
    VAF_ParArray(3)=VAF_par3;
    VAF_ParArray(4)=VAF_par4;

    sFRF_parArray(:,1)=sFRF1_par(1:ftfLim);
    sFRF_parArray(:,2)=sFRF2_par(1:ftfLim);
    sFRF_parArray(:,3)=sFRF3_par(1:ftfLim);
    sFRF_parArray(:,4)=sFRF4_par(1:ftfLim);

    sFRF_array(:,1)=ydata1;
    sFRF_array(:,2)=ydata2;
    sFRF_array(:,3)=ydata3;
    sFRF_array(:,4)=ydata4;

    cohArray(:,1)=cohydata1;
    cohArray(:,2)=cohydata2;
    cohArray(:,3)=cohydata3;
    cohArray(:,4)=cohydata4;

    ftfArray = ftf1(1:ftfLim);
end