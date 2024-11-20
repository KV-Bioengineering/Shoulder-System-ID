%% TrqCalc - Function for calculating Torque from a Stiffness Frequency Response Function (sFRF)
function TrqEst = TrqCalc(sFRF,inAngleData)

    winAngle=tukeywin(length(inAngleData),0.02).*inAngleData; % windowing angle for fft

    fftAngle2s = fft(winAngle,length(sFRF)*2); % 2-sided Angle in freq domain
    fftAngle = fftAngle2s(1:size(fftAngle2s,1)/2); % 1-sided Angle

    fftTrq = fftAngle.*sFRF; % calculate Trq in freq domain
    fftTrq2s = [fftTrq(1); fftTrq(2:end); flipud(conj(fftTrq(2:end)))]; % convert to 2 sided Trq

    TrqEst = ifft(fftTrq2s,length(inAngleData),'symmetric'); % ifft of freq-domain Trq to get time-domain Trq

end