% PECAR_resampshuffletest_frequency.m
%
% This analysis does not generate new data sets. Instead it uses the
% existing data and runs the analysis on resampled trials (alternative
% hypothesis) compared to shuffling the responses for trials (e.g.
% shuffling the right responses only).
%
% Predictions:
% If the relationship between responses is driving the result then this
% should be lost in the shuffle version and there should be no significant
% difference between P1 and P2. 
% 
% A power analysis can be carried out in the same way as before. Are the
% distributions of the two versions overlapping?
%
% 10/1/2017 JCM update from PECAR_resampshuffletest.m to also do frequency analysis
%
% Power analysis
% From Dugue 2015 it sounds like they took the difference between P1 and P2
% for every delay then FFT of this sequence of 13 points (13 delays tested
% in this experiment, 15 in the 2015 paper) to give 6 frequency amplitudes.
% In the permutation test analysis this will be compared with shuffled
% labels on the delays.

clear all; %close all

validity = 1;   % 1 = valid, 2 = invalid

savestuff = 0;  % Do you want to save figures?
seed = 9192017;
rng(seed)
Nsamp = 100;  % Number of samples to bootstrap
scalehist = Nsamp/4;

datadr = 'C:\Users\Kit Moreland\Dropbox\UW\Research\DividedAttention\PECAR_DugueSenoussi\';
datafile = dir([datadr,'datastruct*.mat']);
load(fullfile(datadr, datafile.name))

switch validity
    case 1 % valid
        data = valid;
    case 2 % invalid
        data = invalid;
end

%% Resample and analyse, then shuffle and reanalyse for each run
nsubj = length(data);
ndelays = length(unique(data{1}(:,5)));

P1_resamp = zeros(Nsamp,nsubj,ndelays);
P2_resamp = zeros(Nsamp,nsubj,ndelays);

P1_shuff = zeros(Nsamp,nsubj,ndelays);
P2_shuff = zeros(Nsamp,nsubj,ndelays);

for ss = 1:nsubj
    
    curd = data{ss}; % Pull current subj's data from cell array
    delays = curd(:,5); % delays for this subject
    
    for i = 1:Nsamp
        
        for dl = 1:ndelays
            curdelaydata = find(delays == dl);
            nt = length(curdelaydata);  % Number of trials
            r1 = curd(curdelaydata,3);     % Response 1
            r2 = curd(curdelaydata,4);     % Response 2
            
            %--- Resampling
            % Sample with replacement a sample of the same size for this
            % subject. Need to sample the index so the pairing of probe 1 and
            % two always remains matched.
            idx = datasample(1:nt,nt);
            r1tmp_resamp = r1(idx);
            r2tmp_resamp = r2(idx);
            delaystmp_resamp = delays(idx);
            
            % Do the normal analysis for this sample and store
            [P1_resamp(i,ss,dl), P2_resamp(i,ss,dl)] = quadsolve(r1tmp_resamp, r2tmp_resamp, nt);
            
            %--- Shuffle
            idx_sh = Shuffle(idx);
            r1tmp_sh = r1(idx_sh); % Only shuffle one set of responses
            r2tmp_sh = r2(idx);
            
            % Sanity check: If this shuffling works then it should keep things paired and using the same index should return the mean of the data under the original method.
            %         r1tmp_sh = r1(idx);
            
            % Do the normal analysis for this sample and store
            [P1_shuff(i,ss,dl), P2_shuff(i,ss,dl)] = quadsolve(r1tmp_sh, r2tmp_sh, nt);
        end
    end
end

P1_resamp_mean = squeeze(mean(P1_resamp,1));
P2_resamp_mean = squeeze(mean(P2_resamp,1));
P1_shuff_mean = squeeze(mean(P1_shuff,1));
P2_shuff_mean = squeeze(mean(P2_shuff,1));

Pdiff_resamp_mean = squeeze(mean(P1_resamp - P2_resamp,1));
Pdiff_shuff_mean = squeeze(mean(P1_shuff - P2_shuff,1));

P1_resamp_std = squeeze(std(P1_resamp,1));
P2_resamp_std = squeeze(std(P2_resamp,1));
P1_shuff_std = squeeze(std(P1_shuff,1));
P2_shuff_std = squeeze(std(P2_shuff,1));

P1P2_resamp_diff = P1_resamp_mean - P2_resamp_mean;
P1P2_shuff_diff = P1_shuff_mean - P2_shuff_mean;

% Figure of P1/P2 at different delays
figure(1)
clf
subplot(1,2,1)
hold on
plot(1:ndelays,mean(P1_resamp_mean),'ro-')
plot(1:ndelays,mean(P2_resamp_mean),'bo-')
set(gca,'XTick',1:2:14,'XTickLabel',linspace(40,520,7))
xlabel('Delay (ms)')
ylabel('Probability correct report')
title('Resampled')

subplot(1,2,2)
hold on
plot(1:ndelays,mean(P1_shuff_mean),'ro-')
plot(1:ndelays,mean(P2_shuff_mean),'bo-')
set(gca,'XTick',1:2:14,'XTickLabel',linspace(40,520,7))
xlabel('Delay (ms)')
ylabel('Probability correct report')
title('Permutation')

%% Frequency analysis

t = linspace(40,520,13)/1000; % Need to convert to seconds.

for ss = 1:nsubj
    
    y1 = Pdiff_resamp_mean(ss,:); % time series of current subject
    y2 = Pdiff_shuff_mean(ss,:); % time series of current subject
    
    F1 = fft(y1);
    F2 = fft(y2);
    
    Y1 = complex2real(F1,t);
    Y2 = complex2real(F2,t);
    
    ampresamp(ss,:) = Y1.amp;
    ampshuff(ss,:) = Y2.amp;
    
end

ampresamp_mean = mean(ampresamp,1);
ampshuff_mean = mean(ampshuff,1);

figure(2)
clf

subplot(2,2,1)
plot(t,y1,'b-');
xlabel('Time (s)');
ylabel('Amplitude');
% ylim([0,1])
title('Resampled: P1P2 difference with delay')

subplot(2,2,2)
plot(Y1.freq,Y1.amp,'o-');
set(gca,'XLim',[0,max(Y1.freq)],'XTick',0:2:max(Y1.freq));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Resampled: Amplitude of frequencies')
ylim([0,.3])

subplot(2,2,3)
plot(t,y2,'b-');
xlabel('Time (s)');
ylabel('Amplitude');
% ylim([0,1])
title('Permutation: P1P2 difference with delay')

subplot(2,2,4)
plot(Y2.freq,Y2.amp,'o-');
set(gca,'XLim',[0,max(Y1.freq)],'XTick',0:2:max(Y1.freq));
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('Permutation: Amplitude of frequencies')
ylim([0,.3])



