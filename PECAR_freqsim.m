% PECAR_freqsim.m
%
% Create some fake data which has or does not have a periodicity to it and
% attempt to reconstruct it.
clear all

%% The frequency analysis occurs on the pdif so start by creating a pdif
% vector.
dt_sampled = .04;
maxt_sampled = .520;
t_sampled = dt_sampled:dt_sampled:maxt_sampled;

dt = .01;
maxt = 1;

t = 0:dt:maxt;
nt = length(t);

% Find the locations of the sampled time points in this experiment
t_loc = dt_sampled/dt+1:dt_sampled/dt:find(t==maxt_sampled);
nt_sampled = length(t_loc);

Fs = 1/dt_sampled; % sampling frequency
T = 1/Fs; % period
L = nt_sampled-1; % number of sampled periods
f = Fs*(0:(L/2))/L; 
freqs = f(2:end); % sampled frequencies

nCycles = 5; % Make this the period of the effect
amp = 1;
phase = 30;

pdif = amp*sin(2*pi*t*nCycles-pi*phase/180);

y = fft(pdif);
id = 2:11;
figure(2)
stem(2*abs(y(id))/nt,'fill')

figure(1)
clf
subplot(3,1,1)
hold on
plot(t, pdif)
xlabel('Time (s)');
set(gca,'XTick',0:.1:maxt);
set(gca,'YTick',-1.5:1.5);
ylabel('P1/P2 Difference');
set(gca,'YLim',amp*3*[-1,1]);
%% Now create several runs
nsubj = 15;
sd = 1; % 0 is no noise

for ss = 1:nsubj
    % added noise
    pdif_group(ss,:) = pdif(t_loc) + sd*randn(1,nt_sampled);
    
    %% Run the same analysis as on the real data on this noisy pdif
    F = fft(pdif_group(ss,:)); % Only use the time points sampled in the experiment
    
    A = complex2real(F,t(t_loc));
    
    pdif_amps(ss,:) = A.amp(1:end-1);

end

pdiff_mean = mean(pdif_group);
pdiff_se = std(pdif_group,1)/sqrt(nsubj);

pdif_amps_mean = mean(pdif_amps);
pdif_amps_se = std(pdif_amps)/sqrt(nsubj);

% Plot the mean data of pdif with periodicity
subplot(3,1,2)
hold on
errorbar(t(t_loc), pdiff_mean, pdiff_se, 'ko-')
xlabel('Time (s)');
set(gca,'XTick',0:.1:maxt);
set(gca,'YTick',-1.5:1.5);
ylabel('P1/P2 Difference');
set(gca,'YLim',amp*3*[-1,1]);

% Plot the frequency spectra
subplot(3,1,3)
hold on
errorbar(freqs,pdif_amps_mean,pdif_amps_se,'ko-')
set(gca,'XTick',freqs,'XTickLabel',round(freqs,3,'Significant'))
xlabel('Frequency')
ylabel('Amplitude')
title('Real data: Frequency power')
% Add a line at the real frequency of this fake data
plot([nCycles,nCycles],[0,max(pdif_amps(:))],'k:','LineWidth',2)