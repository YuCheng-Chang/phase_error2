function plv = my_eegPLV(eeg1,eeg2,srate,filtSpec)
%MY_EEGPLV Summary of this function goes here
%   eeg1 and eeg2 are both 1xN
filtPts = fir1(filtSpec.order, 2/srate*filtSpec.range);
eeg1 = filter(filtPts, 1, eeg1, [], 2);
eeg2 = filter(filtPts, 1, eeg2, [], 2);
phase1=angle(hilbert(eeg1));
phase2=angle(hilbert(eeg2));
plv=abs(mean(exp(1i*(phase1-phase2))));
end

