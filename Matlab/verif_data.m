function [] = verif_data(data)
    % This function verify Benford's law on provided data
    linData = squeeze(reshape(data,[1 prod(size(data),'all')]));
    
    fSigNum = fix(abs(linData) .* 10 .^ (-floor(log10(abs(linData)))));
    
    figure();
    histogram(fSigNum);
    hold on
    c=1:0.1:9;
    plot(c,log10(1+1./c).*length(fSigNum),'LineWidth',3);
end