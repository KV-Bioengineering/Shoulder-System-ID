%% TrqConv - Alternate Function for calculating Torque from a Stiffness Frequency Response Function (sFRF)
function TrqEst = TrqConv(sFRF,inAngleData)

    sIRF_nRot = ifft([sFRF;flipud(conj(sFRF(2:end)))]);
    sIRF = ifftshift(sIRF_nRot);

    TrqEst = conv(tukeywin(size(inAngleData,1),0.02).*inAngleData,sIRF,'same');
end