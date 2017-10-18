% PECAR_freqreplication.m

clear all; %close all

validity = 1;   % 1 = valid, 2 = invalid

Fs = 1/.04; % sampling frequency
T = 1/Fs; % period
L = 13-1; % number of sampled periods
f = Fs*(0:(L/2))/L; 
freqs=f(2:end); % sampled frequencies
t = linspace(40,520,13);

seed = 9192017;
rng(seed)

datadr = 'C:\Users\Kit Moreland\Dropbox\UW\Research\DividedAttention\PECAR_DugueSenoussi\PECAR\';
datafile = dir([datadr,'datastruct*.mat']);
load(fullfile(datadr, datafile.name))

switch validity
    case 1 % valid
        data = valid;
        colors = [0,0,1; 0,.5,.8];
    case 2 % invalid
        data = invalid;
        colors = [1,0,0; .6,.4,0];
end

nsubj = length(data);
ndelays = length(unique(data{1}(:,5)));

P1 = zeros(nsubj,ndelays);
P2 = zeros(nsubj,ndelays);
P1_shuff = zeros(nsubj,ndelays);
P2_shuff = zeros(nsubj,ndelays);

for ss = 1:nsubj
    
    curd = data{ss}; % Pull current subj's data from cell array
    delays = curd(:,5); % delays for this subject
    
    for dl = 1:ndelays
        curdelaydata = find(delays == dl);
        nt = length(curdelaydata);  % Number of trials
        r1 = curd(curdelaydata,3);     % Response 1
        r2 = curd(curdelaydata,4);     % Response 2
        
        [P1(ss,dl), P2(ss,dl)] = quadsolve(r1, r2, nt);
        
        
        %--- Shuffle
        r1tmp_sh = Shuffle(r1); % Only shuffle one set of responses
        r2tmp_sh = r2;
        
        % Do the normal analysis
        [P1_shuff(ss,dl), P2_shuff(ss,dl)] = quadsolve(r1tmp_sh, r2tmp_sh, nt);
    end
    
    pdiff(ss,:) = P1(ss,:) - P2(ss,:);
    
    fft_pdiff = pdiff(ss,:);
    A = complex2real(fft_pdiff,t);

    pdif_amps(ss,:) = A.amp(1:end-1);

end

P1_mean = mean(P1,1);
P2_mean = mean(P2,1);
P1_shuff_mean = mean(P1_shuff,1);
P2_shuff_mean = mean(P2_shuff,1);

pdiff_mean = mean(pdiff);
pdiff__shuff_mean = mean(P1_shuff - P2_shuff);

P1_se = std(P1,1)/sqrt(nsubj);
P2_se = std(P2,1)/sqrt(nsubj);
P1_shuff_se = std(P1_shuff,1)/sqrt(nsubj);
P2_shuff_se = std(P2_shuff,1)/sqrt(nsubj);

pdiff_se = std(P1-P2,1)/sqrt(nsubj);
pdiff__shuff_se = std(P1_shuff - P2_shuff)/sqrt(nsubj);

% Figure of Pdiff at different delays
figure(validity)
clf
subplot(1,2,1)
hold on
errorbar(1:ndelays,pdiff_mean,pdiff_se,'ko-')
plot(1:ndelays,zeros(1,ndelays),'--','Color',.25*[1,1,1])
set(gca,'XTick',1:2:14,'XTickLabel',linspace(40,520,7))
% set(gca,'YLim',[0,.8],'YTick',0:.2:.8,'YTickLabel',0:.2:.8)
xlabel('Delay (ms)')
ylabel('P1 - P2')
title('Real data: Pdif')

%% Freq analysis

pdif_amps_mean = mean(pdif_amps);
pdif_amps_se = std(pdif_amps)/sqrt(nsubj);

subplot(1,2,2)
hold on
errorbar(freqs,pdif_amps_mean,pdif_amps_se,'ko-')
set(gca,'XTick',freqs,'XTickLabel',round(freqs,3,'Significant'))
xlabel('Frequency')
ylabel('Amplitude')
title('Real data: Frequency power')

set(gcf,'Position',[2000,700,800,260])

tile
