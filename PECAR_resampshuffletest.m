% PECAR_resampshuffletest.m
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
% 9/19/2017 JC Moreland

clear all; close all

savestuff = 0;  % Do you want to save figures?

seed = 9192017;
rng(seed)
Nsamp = 10000;  % Number of samples to bootstrap
scalehist = Nsamp/2;

datadr = 'C:\Users\Kit Moreland\Dropbox\UW\Research\DividedAttention\PECAR_DugueSenoussi\';
datafile = dir([datadr,'*.mat']);
load(fullfile(datadr, datafile.name))
data = probe_info_valid;

%% Resampling
nsubj = length(data);

P1_resamp = zeros(nsubj,Nsamp);
P2_resamp = zeros(nsubj,Nsamp);

for ss = 1:nsubj
    nt = size(data{ss},1);  % Number of trials
    
    curd = data{ss}; % Pull current subj's data from cell array
    
    r1 = curd(:,3);     % Response 1
    r2 = curd(:,4);     % Response 2
    
    % Resampling and analysis
    for i = 1:Nsamp
        
        % Sample with replacement a sample of the same size for this
        % subject. Need to sample the index so the pairing of probe 1 and
        % two always remains matched.
        idx = datasample(1:nt,nt);
        r1tmp = r1(idx);
        r2tmp = r2(idx);
        
        % Do the normal analysis for this sample and store
        [P1_resamp(ss,i), P2_resamp(ss,i)] = quadsolve(r1tmp, r2tmp, nt);
        
    end
end

P1_resamp_mean = mean(P1_resamp,1);
P2_resamp_mean = mean(P2_resamp,1);

P1P2_resamp_diff = P1_resamp_mean - P2_resamp_mean;

figure(1)
clf
subplot(3,1,1)
histogram(P1P2_resamp_diff,'Binwidth',.01,'FaceColor',[188 143 143]/255)
xlim([0,.5])
xlabel('(P1 - P2)')
ylabel('Frequency')
title('Resampling P1/P2 Difference')
%% Shuffling
nsubj = length(data);

P1_shuff = zeros(nsubj,Nsamp);
P2_shuff = zeros(nsubj,Nsamp);

for ss = 1:nsubj
    nt = size(data{ss},1);  % Number of trials
    
    curd = data{ss}; % Pull current subj's data from cell array
    
    r1 = curd(:,3);     % Response 1
    r2 = curd(:,4);     % Response 2
    
    % Shuffling and analysis
    for i = 1:Nsamp
        
        % Only shuffle one set of responses
        idx = Shuffle(1:nt);
        r1tmp = r1(idx);
        r2tmp = r2;
        
        % Sanity check: If this shuffling works then it should keep things paired and using the same index should return the mean of the data under the original method. 
%         r2tmp = r2(idx);
        
        % Do the normal analysis for this sample and store
        [P1_shuff(ss,i), P2_shuff(ss,i)] = quadsolve(r1tmp, r2tmp, nt);
        
    end
end

P1_shuff_mean = mean(P1_shuff,1);
P2_shuff_mean = mean(P2_shuff,1);

P1P2_shuff_diff = P1_shuff_mean - P2_shuff_mean;

figure(1)
subplot(3,1,2)
histogram(P1P2_shuff_diff,'Binwidth',.01,'FaceColor',[237 145 33]/255)
xlim([0,.5])
xlabel('(P1 - P2)')
ylabel('Frequency')
title('Shuffled P1/P2 Difference')

%% Calculate the achieved power under the two scenarios above

% Plot distributions on the same axis
figure(1)
subplot(3,1,3)
hold on
histogram(P1P2_resamp_diff,'Binwidth',.01,'FaceColor',[188 143 143]/255)
histogram(P1P2_shuff_diff,'Binwidth',.01,'FaceColor',[237 145 33]/255)
xlim([0,.5])
xlabel('(P1 - P2) value')
ylabel('Frequency')
title('Comparison P1/P2 Difference')


p = .05; % What is our alpha?

% Where is the critical value between the two P1-P2 distributions?
crit = prctile(P1P2_shuff_diff,(1-p)*100);

plot([crit,crit],get(gca,'Ylim'),'r-','LineWidth',2)
legend({'Resampled','Shuffled','p-crit'})

power_achieved = sum(P1P2_resamp_diff > crit)/length(P1P2_resamp_diff);
fprintf('Power acheived: %.1f%%\n',power_achieved*100)
ylim([0,scalehist])