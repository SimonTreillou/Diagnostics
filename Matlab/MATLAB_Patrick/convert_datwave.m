load freq_amp_SWASH.txt
data(:,1)=freq_amp_SWASH(:,2);    % amp
data(:,2)=1./freq_amp_SWASH(:,1); % per
data(:,3)=freq_amp_SWASH(:,4);    % pha
save('datwaves.txt','data','-ascii')


