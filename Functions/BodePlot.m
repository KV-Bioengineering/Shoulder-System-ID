%% Graphing
function [bdplt] = BodePlot (sFRF,f,coh,titletxt,sFRF_par) % Parametric FRF vs non-parametric FRF

    bdfig = figure('Name',"Bodeplot - "+titletxt);
    bdplt = tiledlayout("vertical");
    magplt = nexttile;
    phaseplt = nexttile;
    cohplt = nexttile;
    if ~exist("sFRF_par","var")
        % Magnitude v freq
        semilogy(magplt,f,abs(sFRF/size(sFRF,1)))
        % Phase v freq
        plot(phaseplt,f,rad2deg(unwrap(angle(sFRF))))
        % Coh v freq
        semilogx(cohplt,f,coh);
    else
        semilogy(magplt,f,abs(sFRF/size(sFRF,1)),f,abs(sFRF_par/size(sFRF_par,1)))
        % Phase v freq
        plot(phaseplt,f,rad2deg(unwrap(angle(sFRF))),f,rad2deg(unwrap(angle(sFRF_par))))
        % Coh v freq
        plot(cohplt,f,coh);
        legend(magplt,["Non-parametric","Parametric"],'Location','best')
    end
    title(bdplt,titletxt)
    xlabel(bdplt,"Frequency (Hz)")
    ylabel(magplt,"Magnitude (log)")
    ylabel(phaseplt,"Phase (deg)")
    ylabel(cohplt,"Coh")
    % exportgraphics(bplot,"Plots.pdf","Append",true)
end
