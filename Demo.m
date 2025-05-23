%% Ambient noise data
clear

% load an example ambient noise data, including 32 channels
load('Demo_data_ambient.mat')
% In which, 'Uz' is the vertical displacement noise field, 'Time' is the
% recorded time coordinate, 'rloc' is the coordinate of the receiver on the
% xoy plane.

% Perform cross-correlation calculations
% A sliding window with a length of 2000 sampling points and an overlap
% length of 1900 sampling points is used; The results of different
% normalization treatments are calculated.
CC1=mCCFJ.correlate(Uz,rloc,Fs,"AX","xyz","CL",2000,"OL",1900); % No 
CC2=mCCFJ.correlate(Uz,rloc,Fs,"AX","xyz","CL",2000,"OL",1900,"NT","OneBit");
CC3=mCCFJ.correlate(Uz,rloc,Fs,"AX","xyz","CL",2000,"OL",1900,"NF","PSD");
CC4=mCCFJ.correlate(Uz,rloc,Fs,"AX","xyz","CL",2000,"OL",1900,"NF",'ABS');

% Perform wFH transform
% ["Fun","H2","Win","Hamming_1"] is the default configuration of 'mCCFJ.transform'
FJ1=mCCFJ.transform(CC1, 150:500, [1 80]);
FJ2=mCCFJ.transform(CC2, 150:500, [1 80]);
FJ3=mCCFJ.transform(CC3, 150:500, [1 80]);
FJ4=mCCFJ.transform(CC4, 150:500, [1 80]);

% Plot the real part of its dispersion spectrum
figure;
subplot(2,2,1);FJ=FJ1;
imagesc(FJ.frq,FJ.vel,real(FJ.dsp)./max(abs(real(FJ.dsp)),[],1));
set(gca,'YDir','normal');clim([0 1]);
xlabel('Frequency (Hz)');ylabel('Phase velocity (m/s)')
title('Raw data')
hold on
plot(Dispersion_curve.freq, Dispersion_curve.velocity,'-');

subplot(2,2,2);FJ=FJ2;
imagesc(FJ.frq,FJ.vel,real(FJ.dsp)./max(abs(real(FJ.dsp)),[],1));
set(gca,'YDir','normal');clim([0 1]);
xlabel('Frequency (Hz)');ylabel('Phase velocity (m/s)')
title('Onebit time-domain normalization')
hold on
plot(Dispersion_curve.freq, Dispersion_curve.velocity,'-');

subplot(2,2,3);FJ=FJ3;
imagesc(FJ.frq,FJ.vel,real(FJ.dsp)./max(abs(real(FJ.dsp)),[],1));
set(gca,'YDir','normal');clim([0 1]);
xlabel('Frequency (Hz)');ylabel('Phase velocity (m/s)')
title('PSD frequency-domain normalization')
hold on
plot(Dispersion_curve.freq, Dispersion_curve.velocity,'-');

subplot(2,2,4);FJ=FJ4;
imagesc(FJ.frq,FJ.vel,real(FJ.dsp)./max(abs(real(FJ.dsp)),[],1));
set(gca,'YDir','normal');clim([0 1]);
xlabel('Frequency (Hz)');ylabel('Phase velocity (m/s)')
title('ABS frequency-domain normalization')
hold on
plot(Dispersion_curve.freq, Dispersion_curve.velocity,'-');




%% Active source data
clear

% load an example active source data, including 20 channels
load('Demo_data_active.mat')
% In which, 'Uz' is the vertical displacement field, 'Time' is the recorded
% time coordinate, 'rloc' is the coordinate of the receiver on the xoy
% plane.

% For active source data, we do not need to do cross-correlation between
% two station. We can simply bring them directly into the 'mCCFJ.transform'
% function to analyze their dispersion. However, the input object required
% by 'mCCFJ.transform' function is the output of 'mCCFJ.correlate'
% function. The output (named CC) of 'mCCFJ.correlate' at least includes
% frequency domain cross-correlation signal, frequency sequence, and
% inter-station distance. To do this, we build a CC in this format.
CC.ccf=fft(Uz);
CC.ccf=CC.ccf(1:1:ceil(length(Uz(:,1))/2),:); % frequency domain active source data
CC.ccr=sqrt(sum((rloc-[0;0]).^2,1))';  % distance from source, CC.ccr should be a column vector
CC.freq=0:1/Time(end):1/Time(2)/2;     % frequency sequence

% Perform wFH transform
% ["Fun","H2","Win","Hamming_1"] is the default configuration of 'mCCFJ.transform'
FJ=mCCFJ.transform(CC, 150:400, [1 80]); 

% Plot the absolute value of its dispersion spectrum
figure;
imagesc(FJ.frq,FJ.vel,abs(FJ.dsp)./max(abs(FJ.dsp),[],1));
set(gca,'YDir','normal');clim([0 1]);
xlabel('Frequency (Hz)');ylabel('Phase velocity (m/s)')
hold on
plot(Dispersion_curve.freq, Dispersion_curve.velocity,'-');

% Another way to construct CC is to deconvolute all channels to the first
% channel to remove the influence of source time function.
CC1.ccf=CC.ccf(:,2:end)./CC.ccf(:,1); % this deconvolution operation is only illustrative
CC1.ccr=sqrt(sum((rloc(:,2:end)-rloc(:,1)).^2,1))'; % distance from the first channel
CC1.freq=0:1/Time(end):1/Time(2)/2;

% Perform wFH transform
FJ1=mCCFJ.transform(CC1, 150:400, [1 80]);

% Plot the absolute value of its dispersion spectrum
figure;
imagesc(FJ1.frq,FJ1.vel,abs(FJ1.dsp)./max(abs(FJ1.dsp),[],1));
set(gca,'YDir','normal');clim([0 1]);
xlabel('Frequency (Hz)');ylabel('Phase velocity (m/s)')
hold on
plot(Dispersion_curve.freq, Dispersion_curve.velocity,'-');

% If the spatial window is not used, the results are as follows
FJ2=mCCFJ.transform(CC1, 150:400, [1 80],"Win","No");
figure;
imagesc(FJ2.frq,FJ2.vel,abs(FJ2.dsp)./max(abs(FJ2.dsp),[],1));
set(gca,'YDir','normal');clim([0 1]);
xlabel('Frequency (Hz)');ylabel('Phase velocity (m/s)')
hold on
plot(Dispersion_curve.freq, Dispersion_curve.velocity,'-');

