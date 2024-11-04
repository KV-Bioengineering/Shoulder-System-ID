function TrqDEst = DeviceTrqCalc(sFRF_D,inAngleData)

    winAngle=tukeywin(length(inAngleData),0.015).*inAngleData; % windowing angle for fft

    % fAngle = fft(inAngleD); % testing
    % fAngleD = fAngle(1:size(fAngle,1)/2+1); % 1-sided Angle - testing

    fftAngle2s = fft(winAngle); % 2-sided Angle in freq domain
    fftAngle = fftAngle2s(1:size(fftAngle2s,1)/2+1); % 1-sided Angle

    fftTrqD = fftAngle.*sFRF_D; % calculate Trq in freq domain
    fftTrq2sD = [fftTrqD(1); fftTrqD(2:end); flipud(conj(fftTrqD(2:end)))]; % convert to 2 sided Trq

    TrqDEst = ifft(fftTrq2sD); % ifft of freq-domain Trq to get time-domain Trq

end