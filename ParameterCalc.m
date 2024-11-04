%% Parameter Calculation
function [Pararray,sFRF_array,sFRF_parArray,cohArray,VAF_NParray,VAF_ParArray,ftfArray] = ParameterCalc(Fs,inAnglData,outTrqData,deviceData)
    Ts = 1/Fs;
    % NSig = 40001; % <- change to signal length in samples
    % fRes = Fs/NSig;
    % f = 0:fRes:500;
    t=(0:10000)*Ts; %time vector
    t=t';
    
    NSplt = 10001; % size of segments
    % nfft = 2^nextpow2(Nwin);
    fRes_Splt = Fs/NSplt;
    f_Splt = 0:fRes_Splt:500;

    % fLim = floor(10/fRes);
    % w = 2*pi*f(1:fLim);
    
    fLim_Splt=floor((10/fRes_Splt)); %Should be up to 10 Hz
    w_Splt=2*pi*f_Splt(1:fLim_Splt);

    %Butterworth Filter initialization
    fc=10; %Cuttoff Freq
    [b,a] = butter(2,fc/(Fs/2),'low');

    %   Device RF -----------------
    Dwin = 3000;
    Dov = 0.6;

    inAngleD = detrend(deviceData(:,1),1);
    outTrqD = detrend(deviceData(:,2),1);

    sFRF_D = tfestimate(inAngleD,outTrqD,Dwin,Dov*Dwin,size(inAngleD,1),Fs);
    cohD = mscohere(inAngleD,outTrqD,Dwin,Dov*Dwin,size(inAngleD,1),Fs);

    %   Angle and Torque ---------------
    inAngl = detrend(inAnglData,1);

    OutTrqDEst = DeviceTrqCalc(sFRF_D,inAngl);

    outTrq = detrend(outTrqData,1) - OutTrqDEst;

    %Split into 4 segments of 10s each

    inAngl1 = inAngl(1:10001);
    inAngl2 = inAngl(10001:20001);
    inAngl3 = inAngl(20001:30001);
    inAngl4 = inAngl(30001:40001);

    outTrq1 = outTrq(1:10001);
    outTrq2 = outTrq(10001:20001);
    outTrq3 = outTrq(20001:30001);
    outTrq4 = outTrq(30001:40001);

    % outTrq1 = detrend(outTrqData(1:10001),0) - DeviceTrqCalc(sFRF_D(1:5001),inAngl1);
    % outTrq2 = detrend(outTrqData(10001:20001),0) - DeviceTrqCalc(sFRF_D(1:5001),inAngl2);
    % outTrq3 = detrend(outTrqData(20001:30001),0) - DeviceTrqCalc(sFRF_D(1:5001),inAngl3);
    % outTrq4 = detrend(outTrqData(30001:40001),0) - DeviceTrqCalc(sFRF_D(1:5001),inAngl4);

    % outTrq1 = detrend(outTrqData(1:10001),1) - DeviceTrqCalc(sFRF_D,inAngl1);
    % outTrq2 = detrend(outTrqData(10001:20001),1) - DeviceTrqCalc(sFRF_D,inAngl2);
    % outTrq3 = detrend(outTrqData(20001:30001),1) - DeviceTrqCalc(sFRF_D,inAngl3);
    % outTrq4 = detrend(outTrqData(30001:40001),1) - DeviceTrqCalc(sFRF_D,inAngl4);


    %   TFestimate ----------------------
    win = []; % window size
    ov = [];
    % ov = 0.8*win; % overlap

    % [sFRF,ftf] = tfestimate(inAngl,outTrq,win,ov.*win,f,Fs);

    [sFRF1,ftf1] = tfestimate(inAngl1,outTrq1,win,ov,f_Splt,Fs); % cross power spectral density between the input and the output
    [sFRF2,ftf2] = tfestimate(inAngl2,outTrq2,win,ov,f_Splt,Fs);
    [sFRF3,ftf3] = tfestimate(inAngl3,outTrq3,win,ov,f_Splt,Fs);
    [sFRF4,ftf4] = tfestimate(inAngl4,outTrq4,win,ov,f_Splt,Fs);

    % % Correct first term of sFRF -> needs to be real number, but there is some minor (10^-15) complex component, probably something to do with number of fft points, or just because real data
    % sFRF1 = [real(sFRF1(1)), sFRF1(2:end)];
    % sFRF2 = [real(sFRF2(1)), sFRF2(2:end)];
    % sFRF3 = [real(sFRF3(1)), sFRF3(2:end)];
    % sFRF4 = [real(sFRF4(1)), sFRF4(2:end)];
    %In Yahya's thesis he wanted to calculate 4 IBK PARS per trial (1 IBK par for each 10s, rejecting the ones that had low VAF

    %   Non-parametric VAF ----------------
    % fftAngle2s = fft(tukeywin(NSig,0.02).*inAngl);

    % fftAngle2s_1 = fft(tukeywin(NSplt,0.02).*inAnglData(1:10001)); % 2-sided Angle in freq domain
    % fftAngle2s_2 = fft(tukeywin(NSplt,0.02).*inAnglData(10001:20001)); % 2-sided Angle in freq domain
    % fftAngle2s_3 = fft(tukeywin(NSplt,0.02).*inAnglData(20001:30001)); % 2-sided Angle in freq domain
    % fftAngle2s_4 = fft(tukeywin(NSplt,0.02).*inAnglData(30001:40001)); % 2-sided Angle in freq domain

    fftAngle2s_1 = fft(tukeywin(NSplt,0.02).*inAngl1); % 2-sided Angle in freq domain
    fftAngle2s_2 = fft(tukeywin(NSplt,0.02).*inAngl2); % 2-sided Angle in freq domain
    fftAngle2s_3 = fft(tukeywin(NSplt,0.02).*inAngl3); % 2-sided Angle in freq domain
    fftAngle2s_4 = fft(tukeywin(NSplt,0.02).*inAngl4); % 2-sided Angle in freq domain

    % fftAngle = fftAngle2s(1:size(fftAngle2s,1)/2+1);

    fftAngle_1 = fftAngle2s_1(1:size(fftAngle2s_1,1)/2+1); % 1-sided Angle
    fftAngle_2 = fftAngle2s_2(1:size(fftAngle2s_2,1)/2+1); % 1-sided Angle
    fftAngle_3 = fftAngle2s_3(1:size(fftAngle2s_3,1)/2+1); % 1-sided Angle
    fftAngle_4 = fftAngle2s_4(1:size(fftAngle2s_4,1)/2+1); % 1-sided Angle

    % fftTrq = fftAngle.*sFRF';
    
    fftTrq_1 = fftAngle_1.*sFRF1'; % calculate Trq in freq domain
    fftTrq_2 = fftAngle_2.*sFRF2';
    fftTrq_3 = fftAngle_3.*sFRF3';
    fftTrq_4 = fftAngle_4.*sFRF4';

    % fftTrq2s = [fftTrq(1); fftTrq(2:end); flipud(conj(fftTrq(2:end)))];

    fftTrq2s_1 = [fftTrq_1(1); fftTrq_1(2:end); flipud(conj(fftTrq_1(2:end)))]; % convert to 2 sided Trq
    fftTrq2s_2 = [fftTrq_2(1); fftTrq_2(2:end); flipud(conj(fftTrq_2(2:end)))];
    fftTrq2s_3 = [fftTrq_3(1); fftTrq_3(2:end); flipud(conj(fftTrq_3(2:end)))];
    fftTrq2s_4 = [fftTrq_4(1); fftTrq_4(2:end); flipud(conj(fftTrq_4(2:end)))];

    % Trq_NP = ifft(fftTrq2s);

    %Take real value of Trq, there seems to be some minor error (order of 10^-15) from sFRF
    Trq_NP1 = real(ifft(fftTrq2s_1)); % ifft of freq-domain Trq to get time-domain Trq
    Trq_NP2 = real(ifft(fftTrq2s_2));
    Trq_NP3 = real(ifft(fftTrq2s_3));
    Trq_NP4 = real(ifft(fftTrq2s_4));

    % Trq_NP1 = ifft(fftTrq2s_1); % ifft of freq-domain Trq to get time-domain Trq
    % Trq_NP2 = ifft(fftTrq2s_2);
    % Trq_NP3 = ifft(fftTrq2s_3);
    % Trq_NP4 = ifft(fftTrq2s_4);

    % figure
    % tiledlayout("flow");
    % nexttile
    % plot(t,outTrq1,t,Trq_NP1);
    % nexttile
    % plot(t,outTrq2,t,Trq_NP2);
    % nexttile
    % plot(t,outTrq3,t,Trq_NP3);
    % nexttile
    % plot(t,outTrq4,t,Trq_NP4);

    % Trq_filt1 = filtfilt(b,a,outTrq(1:10001));
    % Trq_filt2 = filtfilt(b,a,outTrq(10001:20001));
    % Trq_filt3 = filtfilt(b,a,outTrq(20001:30001));
    % Trq_filt4 = filtfilt(b,a,outTrq(30001:40001));

    Trq_filt1 = filtfilt(b,a,outTrq1);
    Trq_filt2 = filtfilt(b,a,outTrq2);
    Trq_filt3 = filtfilt(b,a,outTrq3);
    Trq_filt4 = filtfilt(b,a,outTrq4);

    Trq_NPFilt1 = filtfilt(b,a,Trq_NP1);
    Trq_NPFilt2 = filtfilt(b,a,Trq_NP2);
    Trq_NPFilt3 = filtfilt(b,a,Trq_NP3);
    Trq_NPFilt4 = filtfilt(b,a,Trq_NP4);

    VAF_NP1 = 100*(1-((var(Trq_filt1-Trq_NPFilt1)))./var(Trq_filt1));
    VAF_NP2 = 100*(1-((var(Trq_filt2-Trq_NPFilt2)))./var(Trq_filt2));
    VAF_NP3 = 100*(1-((var(Trq_filt3-Trq_NPFilt3)))./var(Trq_filt3));
    VAF_NP4 = 100*(1-((var(Trq_filt4-Trq_NPFilt4)))./var(Trq_filt4));

    % figure
    % tiledlayout("flow");
    % nexttile
    % plot(t,Trq_filt1,t,Trq_NPFilt1);
    % nexttile
    % plot(t,Trq_filt2,t,Trq_NPFilt2);
    % nexttile
    % plot(t,Trq_filt3,t,Trq_NPFilt3);
    % nexttile
    % plot(t,Trq_filt4,t,Trq_NPFilt4);

    %   Coherence ------------------------
    coh1 = mscohere(inAngl1,outTrq1,win,ov,f_Splt,Fs);
    coh2 = mscohere(inAngl2,outTrq2,win,ov,f_Splt,Fs);
    coh3 = mscohere(inAngl3,outTrq3,win,ov,f_Splt,Fs);
    coh4 = mscohere(inAngl4,outTrq4,win,ov,f_Splt,Fs);

    % coh3 = mscohere(inAngl4,outTrq3,win,ov,f,Fs); % for testing if VAF is less than 50%

    % figure
    % plot(angle(sFRF(1:fLim)))

    % Cross & auto power spectrum testing
    % [Puy,ftf] = cpsd(inAngl,outTrq,[],[],Nwin,Fs); % cross power spectral density between the input and the output
    % Puu = cpsd(inAngl,inAngl,[],[],Nwin,Fs); % auto power spectral density of the input
    % Pyy = cpsd(outTrq,outTrq,[],[],Nwin,Fs); % auto power spectral density of the output
    % Pyu = cpsd(outTrq,inAngl,[],[],Nwin,Fs); % cross power spectral density between the input and the output
    % 
    % sFRF = Puy./Puu;

    % cohf=((abs(Puy)).^2)./(Puu.*Pyy);   % Coherence squared
    % cohY=((abs(Pyu)).^2)./(Puu.*Pyy); % Yahya Coherence
    %All give same answer

    % figure
    % plot(ftf(1:fLim),cohm(1:fLim),ftf(1:fLim),coh(1:fLim),ftf(1:fLim),cohY(1:fLim))
    % plot(ftf(1:fLim),cohY(1:fLim),ftf(1:fLim),coh(1:fLim))
    % legend("Matlab Coh","Formula Coh","Yahya Coh")
    % legend("Yahya Coh","Formula Coh")

    %   Parameter Estimation ---------------

    E0=[0.2,2,20];
    opts = optimoptions("lsqnonlin","Algorithm","levenberg-marquardt","StepTolerance",1.000000e-20,"Display","iter","MaxFunctionEvaluations",6000,"MaxIterations",1000);

    %Y-data for parameter estimation
    ydata1 = sFRF1(1:fLim_Splt);
    ydata2 = sFRF2(1:fLim_Splt);
    ydata3 = sFRF3(1:fLim_Splt);
    ydata4 = sFRF4(1:fLim_Splt);

    cohydata1 = coh1(1:fLim_Splt);
    cohydata2 = coh2(1:fLim_Splt);
    cohydata3 = coh3(1:fLim_Splt);
    cohydata4 = coh4(1:fLim_Splt);

    %Cost Function(s):

    % err = @(E) (abs(ydata)+angle(ydata)) - (abs(E(1)*((i*w).^2)+E(2)*i*w+E(3))+angle(E(1)*((i*w).^2)+E(2)*i*w+E(3))); % textbook fun
    % err = @(E) cohydata.*(abs(ydata)+angle(ydata)) - (abs(E(1)*((i*w).^2)+E(2)*i*w+E(3))+angle(E(1)*((i*w).^2)+E(2)*i*w+E(3))); %coh weighting on textbook fun
    % err = @(E) [cohydata.*(abs(ydata)-abs((E(1)*((i*w).^2)+E(2)*i*w+E(3)))); cohydata.*(angle(ydata)-angle((E(1)*((i*w).^2)+E(2)*i*w+E(3))))];
    % err = @(E) [(abs(ydata)-abs((E(1)*((i*w).^2)+E(2)*i*w+E(3)))); (angle(ydata)-angle((E(1)*((i*w).^2)+E(2)*i*w+E(3))))];

    % err1 = @(E1) [cohydata1.*(abs(ydata1)-abs((E1(1)*((i*w_Splt).^2)+E1(2)*i*w_Splt+E1(3)))); cohydata1.*(angle(ydata1)-angle((E1(1)*((i*w_Splt).^2)+E1(2)*i*w_Splt+E1(3))))];
    % err2 = @(E2) [cohydata2.*(abs(ydata2)-abs((E2(1)*((i*w_Splt).^2)+E2(2)*i*w_Splt+E2(3)))); cohydata2.*(angle(ydata2)-angle((E2(1)*((i*w_Splt).^2)+E2(2)*i*w_Splt+E2(3))))];
    % err3 = @(E3) [cohydata3.*(abs(ydata3)-abs((E3(1)*((i*w_Splt).^2)+E3(2)*i*w_Splt+E3(3)))); cohydata3.*(angle(ydata3)-angle((E3(1)*((i*w_Splt).^2)+E3(2)*i*w_Splt+E3(3))))];
    % err4 = @(E4) [cohydata4.*(abs(ydata4)-abs((E4(1)*((i*w_Splt).^2)+E4(2)*i*w_Splt+E4(3)))); cohydata4.*(angle(ydata4)-angle((E4(1)*((i*w_Splt).^2)+E4(2)*i*w_Splt+E4(3))))];
    
    % err1 = @(E1) [cohydata1.*(abs(ydata1)-abs((E1(1)*((i*w_Splt).^2)+E1(2)*i*w_Splt+E1(3)))); (angle(ydata1)-angle((E1(1)*((i*w_Splt).^2)+E1(2)*i*w_Splt+E1(3))))];
    % err2 = @(E2) [cohydata2.*(abs(ydata2)-abs((E2(1)*((i*w_Splt).^2)+E2(2)*i*w_Splt+E2(3)))); (angle(ydata2)-angle((E2(1)*((i*w_Splt).^2)+E2(2)*i*w_Splt+E2(3))))];
    % err3 = @(E3) [cohydata3.*(abs(ydata3)-abs((E3(1)*((i*w_Splt).^2)+E3(2)*i*w_Splt+E3(3)))); (angle(ydata3)-angle((E3(1)*((i*w_Splt).^2)+E3(2)*i*w_Splt+E3(3))))];
    % err4 = @(E4) [cohydata4.*(abs(ydata4)-abs((E4(1)*((i*w_Splt).^2)+E4(2)*i*w_Splt+E4(3)))); (angle(ydata4)-angle((E4(1)*((i*w_Splt).^2)+E4(2)*i*w_Splt+E4(3))))];

    err1 = @(E1) cohydata1.*(abs(ydata1)+angle(ydata1) - (abs(E1(1)*((i*w_Splt).^2)+E1(2)*i*w_Splt+E1(3))+angle(E1(1)*((i*w_Splt).^2)+E1(2)*i*w_Splt+E1(3))));
    err2 = @(E2) cohydata2.*(abs(ydata2)+angle(ydata2) - (abs(E2(1)*((i*w_Splt).^2)+E2(2)*i*w_Splt+E2(3))+angle(E2(1)*((i*w_Splt).^2)+E2(2)*i*w_Splt+E2(3))));
    err3 = @(E3) cohydata3.*(abs(ydata3)+angle(ydata3) - (abs(E3(1)*((i*w_Splt).^2)+E3(2)*i*w_Splt+E3(3))+angle(E3(1)*((i*w_Splt).^2)+E3(2)*i*w_Splt+E3(3))));
    err4 = @(E4) cohydata4.*(abs(ydata4)+angle(ydata4) - (abs(E4(1)*((i*w_Splt).^2)+E4(2)*i*w_Splt+E4(3))+angle(E4(1)*((i*w_Splt).^2)+E4(2)*i*w_Splt+E4(3))));

    %Least squares fit to estimate IBK
    % E = lsqnonlin(err,E0,[],[],options);
    [E1,resnorm1,residual1] = lsqnonlin(err1,E0,[],[],opts);
    [E2,resnorm2,residual2] = lsqnonlin(err2,E0,[],[],opts);
    [E3,resnorm3,residual3] = lsqnonlin(err3,E0,[],[],opts);
    [E4,resnorm4,residual4] = lsqnonlin(err4,E0,[],[],opts);

    %Calculate angular velocity, damping coefficient, TBD need to pass through to pararray - In Yahya's code this part isn't used in this function and fn and zeta are calculated again in main code
    % I=E(1);
    % B=E(2);
    % K=E(3);
    % Wn=sqrt(K/I);
    % fn=Wn/(2*pi);
    % zeta=B/(2*sqrt(I*K));

%% Parametric Model and VAF - Need to finish stuff below this
    
    %   Parametric functions - in time domain
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

    Trq_sim1=conv(tukeywin(NSplt,0.015).*inAngl1,sIRF1_par,'same'); % Simulated Intrinsic torque
    Trq_sim2=conv(tukeywin(NSplt,0.015).*inAngl2,sIRF2_par,'same');
    Trq_sim3=conv(tukeywin(NSplt,0.015).*inAngl3,sIRF3_par,'same');
    Trq_sim4=conv(tukeywin(NSplt,0.015).*inAngl4,sIRF4_par,'same');

    %testing graphs
    % figure
    % tiledlayout("flow");
    % nexttile
    % plot(tukeywin(NSplt,0.02).*inAngl1);
    % nexttile
    % plot(tukeywin(NSplt,0.02).*inAngl2);
    % nexttile
    % plot(tukeywin(NSplt,0.02).*inAngl3);
    % nexttile
    % plot(tukeywin(NSplt,0.02).*inAngl4);

    % figure
    % tiledlayout("flow");
    % nexttile
    % plot(Trq_sim1);
    % nexttile
    % plot(Trq_sim2);
    % nexttile
    % plot(Trq_sim3);
    % nexttile
    % plot(Trq_sim4);


    % %   Parametric functions - in freq Domain
    % sFRF1_par = E1(1)*((i*w).^2)+E1(2)*i*w+E1(3);
    % sFRF2_par = E2(1)*((i*w).^2)+E2(2)*i*w+E2(3);
    % sFRF3_par = E3(1)*((i*w).^2)+E3(2)*i*w+E3(3);
    % sFRF4_par = E4(1)*((i*w).^2)+E4(2)*i*w+E4(3);
    % 
    % % winAngle=tukeywin(length(inAngleD),0.004).*inAngleD; % windowing angle for fft
    % 
    % fftAngle2s_1 = fft(inAngl1); % 2-sided Angle in freq domain
    % fftAngle2s_2 = fft(inAngl2); % 2-sided Angle in freq domain
    % fftAngle2s_3 = fft(inAngl3); % 2-sided Angle in freq domain
    % fftAngle2s_4 = fft(inAngl4); % 2-sided Angle in freq domain
    % 
    % fftAngle_1 = fftAngle2sD_1(1:size(fftAngle2sD_1,1)/2+1); % 1-sided Angle
    % fftAngle_2 = fftAngle2sD_2(1:size(fftAngle2sD_2,1)/2+1); % 1-sided Angle
    % fftAngle_3 = fftAngle2sD_3(1:size(fftAngle2sD_3,1)/2+1); % 1-sided Angle
    % fftAngle_4 = fftAngle2sD_4(1:size(fftAngle2sD_4,1)/2+1); % 1-sided Angle
    % 
    % fftTrq_1 = fftAngleD_1.*sFRF1_par'; % calculate Trq in freq domain
    % fftTrq_2 = fftAngleD_2.*sFRF2_par';
    % fftTrq_3 = fftAngleD_3.*sFRF3_par';
    % fftTrq_4 = fftAngleD_4.*sFRF4_par';
    % 
    % fftTrq2s_1 = [fftTrqD_1(1); fftTrqD_1(2:end); flipud(conj(fftTrqD_1(2:end)))]; % convert to 2 sided Trq
    % fftTrq2s_2 = [fftTrqD_2(1); fftTrqD_2(2:end); flipud(conj(fftTrqD_2(2:end)))];
    % fftTrq2s_3 = [fftTrqD_3(1); fftTrqD_3(2:end); flipud(conj(fftTrqD_3(2:end)))];
    % fftTrq2s_4 = [fftTrqD_4(1); fftTrqD_4(2:end); flipud(conj(fftTrqD_4(2:end)))];
    % 
    % Trq_sim1 = ifft(fftTrq2sD_1); % ifft of freq-domain Trq to get time-domain Trq
    % Trq_sim2 = ifft(fftTrq2sD_2);
    % Trq_sim3 = ifft(fftTrq2sD_3);
    % Trq_sim4 = ifft(fftTrq2sD_4);

    % Filtered Torques

    % Trq_filt1 = filtfilt(b,a,outTrq1);
    % Trq_filt2 = filtfilt(b,a,outTrq2);
    % Trq_filt3 = filtfilt(b,a,outTrq3);
    % Trq_filt4 = filtfilt(b,a,outTrq4);

    Trq_simFilt1 = filtfilt(b,a,Trq_sim1);
    Trq_simFilt2 = filtfilt(b,a,Trq_sim2);
    Trq_simFilt3 = filtfilt(b,a,Trq_sim3);
    Trq_simFilt4 = filtfilt(b,a,Trq_sim4);

    % figure
    % tiledlayout("flow")
    % nexttile
    % plot(t,Trq_filt1,t,Trq_simFilt1);
    % nexttile
    % plot(t,Trq_filt2,t,Trq_simFilt2);
    % nexttile
    % plot(t,Trq_filt3,t,Trq_simFilt3);
    % nexttile
    % plot(t,Trq_filt4,t,Trq_simFilt4);


    %Parametric model Variance Accounted For

    VAF_par1 = 100*(1-((var(Trq_filt1-Trq_simFilt1)))./var(Trq_filt1));
    VAF_par2 = 100*(1-((var(Trq_filt2-Trq_simFilt2)))./var(Trq_filt2));
    VAF_par3 = 100*(1-((var(Trq_filt3-Trq_simFilt3)))./var(Trq_filt3));
    VAF_par4 = 100*(1-((var(Trq_filt4-Trq_simFilt4)))./var(Trq_filt4));

    %Exclude if VAF_par is less than 50%
    E_tot = [0,0,0];
    nE=0; % number of parameter sets taken into account

    if(VAF_NP1) > 50
        if(VAF_par1) > 85
            E_tot = E_tot+E1;
            nE=nE+1;
        end
    end

    if(VAF_NP2) > 50
        if(VAF_par2) > 85
            E_tot = E_tot+E2;
            nE=nE+1;
        end
    end

    if(VAF_NP3) > 50
        if(VAF_par3) > 85
            E_tot = E_tot+E3;
            nE=nE+1;
        end
    end

    if(VAF_NP4) > 50
        if(VAF_par4) > 85
            E_tot = E_tot+E4;
            nE=nE+1;
        end
    end

E=E_tot/nE;


    % nexttile
    % BodePlot(sFRF1(1:fLim_Splt),ftf1,coh1,"Participant Segment 1",sFRF_par(1:fLim_Splt))
    % nexttile
    % BodePlot(sFRF2_par)
    % nexttile
    % BodePlot(sFRF3_par)
    % nexttile
    % BodePlot(sFRF4_par)
    
    %Parameter Array
    Pararray(1:3)=E;
    % Pararray(1,:)=E1;
    % Pararray(2,:)=E2;
    % Pararray(3,:)=E3;
    % Pararray(4,:)=E4;

    VAF_NParray(1)=VAF_NP1;
    VAF_NParray(2)=VAF_NP2;
    VAF_NParray(3)=VAF_NP3;
    VAF_NParray(4)=VAF_NP4;

    VAF_ParArray(1)=VAF_par1;
    VAF_ParArray(2)=VAF_par2;
    VAF_ParArray(3)=VAF_par3;
    VAF_ParArray(4)=VAF_par4;

    sFRF_array(:,1)=ydata1;
    sFRF_array(:,2)=ydata2;
    sFRF_array(:,3)=ydata3;
    sFRF_array(:,4)=ydata4;

    sFRF_parArray(:,1)=sFRF1_par(1:fLim_Splt);
    sFRF_parArray(:,2)=sFRF2_par(1:fLim_Splt);
    sFRF_parArray(:,3)=sFRF3_par(1:fLim_Splt);
    sFRF_parArray(:,4)=sFRF4_par(1:fLim_Splt);

    cohArray(:,1)=cohydata1;
    cohArray(:,2)=cohydata2;
    cohArray(:,3)=cohydata3;
    cohArray(:,4)=cohydata4;

    ftfArray(:,1)=ftf1(1:fLim_Splt);
    ftfArray(:,2)=ftf2(1:fLim_Splt);
    ftfArray(:,3)=ftf3(1:fLim_Splt);
    ftfArray(:,4)=ftf4(1:fLim_Splt);

end