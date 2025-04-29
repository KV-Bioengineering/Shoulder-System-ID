%% PCalc - function that estimates the parameters of the GH joint, as well as calculates the parametric Variance accounted for
function [Pararray,MnTrq,VAF_NPArray,VAF_ParArray,sFRF_parArray,sFRF_array,cohArray,ftfArray,Trq_filt,Trq_simFilt] = PCalc(Fs,inAnglData,outTrqData,sFRF_D,E0)
    Ts = 1/Fs;
    t=(0:10000)*Ts; %time vector
    t=t';

    NSplt = size(inAnglData,1); % size of segments

    %Butterworth Filter initialization
    fc=10; %Cuttoff Freq
    [b,a] = butter(2,fc/(Fs/2),'low');

    %   Angle & Torque Data -----------------
    % Split Data into 4 data sets
    % Angle Data
    inAngle = inAnglData;

    % Detrended Angle Data
    dtAngle = detrend(inAnglData,1);

    % Device Torque
    dTrq = TrqCalc(sFRF_D,inAngle);

    %Alternate device torque calculation
    % dTrq = TrqConv(sFRF_D,inAngle);

    % Output Torque (Contraction Torque)
    outTrq = outTrqData - dTrq;
    % Detrended Contraction Torque
    dtTrq = detrend(outTrq,1);

    %%  Transfer Function Estimation -----------------
    win = []; % window size
    ov = [];
    % ov = 0.8*win; % for custom overlap

    [sFRF,ftf] = tfestimate(dtAngle,dtTrq,win,ov,NSplt,Fs); % cross power spectral density between the input and the output

    %   Coherence ------------------------
    coh = mscohere(dtAngle,dtTrq,win,ov,NSplt,Fs);

    %% Non-Parametric VAF

    fftAngle2s = fft(tukeywin(NSplt,0.02).*dtAngle); % 2-sided Angle in freq domain
    fftAngle = fftAngle2s(1:size(fftAngle2s,1)/2+1); % 1-sided Angle 

    fftTrq = fftAngle.*sFRF; % calculate Trq in freq domain
    fftTrq2s = [fftTrq(1); fftTrq(2:end); flipud(conj(fftTrq(2:end)))]; % convert to 2 sided Trq

    Trq_NP = ifft(fftTrq2s); % ifft of freq-domain Trq to get time-domain Trq

    Trq_filt = filtfilt(b,a,dtTrq);
    Trq_NPFilt = filtfilt(b,a,Trq_NP);

    VAF_NP = 100*(1-((var(Trq_filt-Trq_NPFilt)))./var(Trq_filt));

    %%   Parameter Estimation ---------------

    if ~exist("E0","var")
        E0=[0.2,2,20];
    end

    opts = optimoptions("lsqnonlin","Algorithm","levenberg-marquardt","StepTolerance",1.000000e-20,"Display","iter","MaxFunctionEvaluations",6000,"MaxIterations",1000);

    %Y-data for parameter estimation
    ftfLim = find(ftf<10,1,'last');

    w_Splt = 2*pi*ftf(1:ftfLim);

    ydata = sFRF(1:ftfLim);

    cohydata = coh(1:ftfLim);

    % err = @(E) cohydata.*(abs(ydata)+angle(ydata) - (abs(E(1)*((i*w_Splt).^2)+E(2)*i*w_Splt+E(3))+angle(E(1)*((i*w_Splt).^2)+E(2)*i*w_Splt+E(3))));

    err = @(E) [cohydata.*(abs(ydata)-abs((E(1)*((i*w_Splt).^2)+E(2)*i*w_Splt+E(3)))); cohydata.*(angle(ydata)-angle((E(1)*((i*w_Splt).^2)+E(2)*i*w_Splt+E(3))))];

    %Least squares fit to estimate IBK
    % E = lsqnonlin(err,E0,[],[],options);
    [E,resnorm1,residual] = lsqnonlin(err,E0,[],[],opts);

    %% Parametric Model and VAF

    cIRF_par = ((2*(-i)*sin(i*(t*(E(2)^2 - 4*E(1)*E(3))^(1/2))/(2*E(1))).*exp(-(E(2)*t)/(2*E(1))))/(E(2)^2 - 4*E(1)*E(3))^(1/2));  %Simulated Compliance IRF (Parametric) <- from the inverse Laplace of the compliance formula

    cFRF_par=fft(cIRF_par)/Fs;   % Parametric Compliance FRF
    sFRF_par=1./cFRF_par;        % Parametric Stiffness FRF
    sIRF_parNrot=ifft(sFRF_par);
    sIRF_par=ifftshift(sIRF_parNrot);

    Trq_sim=conv(tukeywin(NSplt,0.02).*dtAngle,sIRF_par,'same'); % Simulated Intrinsic torque

    %Filtered Simulated Torque
    Trq_simFilt = filtfilt(b,a,Trq_sim);

    %Parametric model Variance Accounted For
    VAF_par = 100*(1-((var(Trq_filt-Trq_simFilt)))./var(Trq_filt));

    %% Function output
    % Pararray contains IBK values of each set
    Pararray = E';
    MnTrq = mean(outTrq);
    VAF_NPArray=VAF_NP;
    VAF_ParArray=VAF_par;
    sFRF_parArray=sFRF_par(1:ftfLim);
    sFRF_array=ydata;
    cohArray=cohydata;
    ftfArray = ftf(1:ftfLim);
end