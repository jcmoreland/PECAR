% PECAR_bootsim2.m
%
% This simulation is an update from PECAR_bootsim to generate groups of
% observers that are different (rather than using the experiment's
% estimates)
%
% 9/12/2017 JC Moreland

clear all; close all

savestuff = 0;  % Do you want to save figures?

seed = 9102017;
rng(seed)
Nboot = 10000;  % Number of samples to bootstrap
Nsim = 100;     % Number of times to simulate the data set
nTrials = 1400; % Nunber of trials per simulation
scalehist = Nsim/2;

datadr = 'C:\Users\Kit Moreland\Dropbox\UW\Research\DividedAttention\PECAR_DugueSenoussi\';
datafile = dir([datadr,'*.mat']);
load(fullfile(datadr, datafile.name))
data = probe_info_valid;

%% Bootstrap P1 and P2 for each observer
nsubj = length(data);

P1 = zeros(1,nsubj);
P1se = zeros(1,nsubj);
P2 = zeros(1,nsubj);
P2se = zeros(1,nsubj);

pc = zeros(nsubj,2);
pcse = zeros(nsubj,2);

for ss = 1:nsubj
    nt = size(data{ss},1);  % Number of trials
    
    curd = data{ss}; % Pull current subj's data from cell array
    
    r1 = curd(:,3);
    r2 = curd(:,4);
    
    % Bootstrap
    P1tmp = zeros(1,Nboot);
    P2tmp = zeros(1,Nboot);
    pctmp = zeros(Nboot,2);
    for bt = 1:Nboot
        % Sample with replacement a sample of the same size for this
        % subject. Need to sample the index so the pairing of probe 1 and
        % two always remains matched.
        idx = datasample(1:nt,nt);
        r1tmp = r1(idx);
        r2tmp = r2(idx);
        
        % Do the normal analysis for this sample and store
        [P1tmp(bt), P2tmp(bt)] = quadsolve(r1tmp, r2tmp, nt);
        % Percent correct fot this sim
        pctmp(bt,:) = mean([r1tmp,r2tmp]);
        
    end
    
    P1(ss) = mean(P1tmp);
    P1se(ss) = std(P1tmp);
    P2(ss) = mean(P2tmp);
    P2se(ss) = std(P2tmp);
    
    pc(ss,:) = mean(pctmp);
    pcse(ss,:) = std(pctmp);
end

% Group P1/P2 means and SEs
P1group = mean(P1);
P1segroup = std(P1)/sqrt(nsubj);
P2group = mean(P2);
P2segroup = std(P2)/sqrt(nsubj);

% Within subject mean P1/P2 SE and variability
P1sewithin = mean(P1se);
P2sewithin = mean(P2se);

% Group PC means and SEs
pcgroup = mean(pc);
pcsegroup = std(pc)/sqrt(nsubj);
% Within subject mean PC SE and variability
pcsewithin = mean(pcse);

figure(1)
clf
subplot(4,3,1)
hold on
plot(repmat([1,2],nsubj,1)', [P1;P2], 'o-', 'Color', .5*[1,1,1])
errorbar([1,2], [mean(P1),mean(P2)], ...
    [std(P1)/sqrt(nsubj),std(P2)/sqrt(nsubj)], ...
    [std(P1)/sqrt(nsubj),std(P2)/sqrt(nsubj)], 'k','Linewidth',2)
xlim([.5,2.5])
ylim([0,1])
set(gca,'XTick',[1,2],'XTickLabel',{'P1','P2'},'Linewidth',1.5,'FontSize',14)
ylabel('Probability correct report')
title('Real data P1 and P2')

% P1
subplot(4,3,4)
histogram(P1,'BinWidth',.05)
xlabel('P1 value')
ylabel('Frequency')
xlim([0,1])
ylim([0,nsubj])
% P2
subplot(4,3,7)
histogram(P2,'BinWidth',.05)
xlabel('P2 value')
ylabel('Frequency')
xlim([0,1])
ylim([0,nsubj])

% pDif
subplot(4,3,10)
histogram((P1-P2),'BinWidth',.05)
xlabel('(P1 - P2) value')
ylabel('Frequency')
xlim([0,1])
ylim([0,nsubj])

%% Use the group variability and mean P1 and P2 along with the within 
% subject variability to simulate a data set under the two hypotheses.

%% 1. Switching
% In every simulation a P1 and P2 is drawn for each subject and the
% analysis run as in the original.

P1_sim_switch = nan(1,Nsim);
P1se_sim_switch = nan(1,Nsim);
P2_sim_switch = nan(1,Nsim);
P2se_sim_switch = nan(1,Nsim);

pc_sim_switch = nan(Nsim,2);
pcse_sim_switch = nan(Nsim,2);

for n = 1:Nsim

    P1tmp = zeros(1,nsubj);
    P2tmp = zeros(1,nsubj);
    pctmp = zeros(nsubj,2);
    
    for ss = 1:nsubj
        % For each subject draw a value of P1 and P2 from the population
        % P1/P2 distribution estimated with the se
        p1_curmean = P1segroup*randn(1) + P1group;
        p2_curmean = P2segroup*randn(1) + P2group;
        
        resp = nan(nTrials,2);
        for i = 1:nTrials
            
            % To get the within subject error resample the P1/P2 using
            % within error.
            p1_cur = P1sewithin*randn(1) + p1_curmean;
            p2_cur = P2sewithin*randn(1) + p2_curmean;

            resp(i,1) = binornd(1,p1_cur);
            resp(i,2) = binornd(1,p2_cur);
            
            resp(i,:) = Shuffle(resp(i,:)); % Shuffle so that these responses could have been on either side
        end
        
        [P1tmp(ss), P2tmp(ss)] = quadsolve(resp(:,1), resp(:,2), nTrials);
        
        pctmp(ss,:) = mean(resp);
    end
    
    % Get the mean P1 and P2 across the simulated observers for this simulation run
    P1_sim_switch(n) = mean(P1tmp);
    P1se_sim_switch(n) = std(P1tmp)/sqrt(nsubj);
    P2_sim_switch(n) = mean(P2tmp);
    P2se_sim_switch(n) = std(P2tmp)/sqrt(nsubj);
    
    pc_sim_switch(n,:) = mean(pctmp);
    pcse_sim_switch(n,:) = std(pctmp)/sqrt(nsubj); % This is the se of the subjects performance
    
end

figure(1)
subplot(4,3,2)
hold on
% Percent correct
errorbar([1,2], mean(pc_sim_switch), ...
    [mean(pcse_sim_switch)], ...
    [mean(pcse_sim_switch)], 'rx','Linewidth',2)
% P1 and P2
plot(repmat([1,2],n,1)', [P1_sim_switch;P2_sim_switch], 'o-', 'Color', .5*[1,1,1])
errorbar([1,2], [mean(P1_sim_switch),mean(P2_sim_switch)], ...
    [mean(P1se_sim_switch),mean(P2se_sim_switch)], ...
    [mean(P1se_sim_switch),mean(P2se_sim_switch)], 'ko','Linewidth',2)

xlim([.5,2.5])
ylim([0,1])
set(gca,'XTick',[1,2],'XTickLabel',{'P1','P2'},'Linewidth',1.5,'FontSize',14)
ylabel('Probability correct report')
title('Simulated data P1 and P2: Switching')

% P1
subplot(4,3,5)
histogram(P1_sim_switch,'BinWidth',.01)
xlabel('P1 value')
ylabel('Frequency')
xlim([0,1])
ylim([0,scalehist])
% P2
subplot(4,3,8)
histogram(P2_sim_switch,'BinWidth',.01)
xlabel('P2 value')
ylabel('Frequency')
xlim([0,1])
ylim([0,scalehist])

% pDif
Pdif_switch = (P1_sim_switch-P2_sim_switch);
subplot(4,3,11)
histogram(Pdif_switch,'BinWidth',.01)
xlabel('(P1 - P2) value')
ylabel('Frequency')
xlim([0,1])
ylim([0,scalehist])

%% 2. Unlimited parallel

P1_sim_unlim = nan(1,Nsim);
P1se_sim_unlim = nan(1,Nsim);
P2_sim_unlim = nan(1,Nsim);
P2se_sim_unlim = nan(1,Nsim);

pc_sim_unlim = nan(Nsim,2);
pcse_sim_unlim = nan(Nsim,2);

for n = 1:Nsim

    P1tmp = zeros(1,nsubj);
    P2tmp = zeros(1,nsubj);
    pctmp = zeros(nsubj,2);
    
    for ss = 1:nsubj
            
        % For each subject draw a percent correct value
        p1_curmean = pcsegroup(1)*randn(1) + pcgroup(1);
        p2_curmean = pcsegroup(2)*randn(1) + pcgroup(2);
        
        resp = nan(nTrials,2);
        for i = 1:nTrials
            
            % To get the within subject error resample the pc using
            % within error.
            p1_cur = pcsewithin(1)*randn(1) + p1_curmean;
            p2_cur = pcsewithin(1)*randn(1) + p2_curmean;
            
            resp(i,1) = binornd(1,p1_cur);
            resp(i,2) = binornd(1,p2_cur);
            
            resp(i,:) = Shuffle(resp(i,:)); % Shuffle so that these responses could have been on either side
        end
        
        [P1tmp(ss), P2tmp(ss)] = quadsolve(resp(:,1), resp(:,2), nTrials);
        
        pctmp(ss,:) = mean(resp);
    end
    
    % Get the mean P1 and P2 across the simulated observers for this simulation run
    P1_sim_unlim(n) = mean(P1tmp);
    P1se_sim_unlim(n) = std(P1tmp)/sqrt(nsubj);
    P2_sim_unlim(n) = mean(P2tmp);
    P2se_sim_unlim(n) = std(P2tmp)/sqrt(nsubj);
    
    pc_sim_unlim(n,:) = mean(pctmp);
    pcse_sim_unlim(n,:) = std(pctmp)/sqrt(nsubj); % This is the se of the subjects performance
    
end


figure(1)
subplot(4,3,3)
hold on
% Percent correct
errorbar([1,2], mean(pc_sim_unlim), ...
    [mean(pcse_sim_unlim)], ...
    [mean(pcse_sim_unlim)], 'rx','Linewidth',2)
% P1 and P2
plot(repmat([1,2],n,1)', [P1_sim_unlim;P2_sim_unlim], 'o-', 'Color', .5*[1,1,1])
errorbar([1,2], [mean(P1_sim_unlim),mean(P2_sim_unlim)], ...
    [mean(P1se_sim_unlim),mean(P2se_sim_unlim)], ...
    [mean(P1se_sim_unlim),mean(P2se_sim_unlim)], 'ko','Linewidth',2)

xlim([.5,2.5])
ylim([0,1])
set(gca,'XTick',[1,2],'XTickLabel',{'P1','P2'},'Linewidth',1.5,'FontSize',14)
ylabel('Probability correct report')
title('Simulated data P1 and P2: Unlimited')

% P1
subplot(4,3,6)
histogram(P1_sim_unlim,'BinWidth',.01)
xlabel('P1 value')
ylabel('Frequency')
xlim([0,1])
ylim([0,scalehist])
% P2
subplot(4,3,9)
histogram(P2_sim_unlim,'BinWidth',.01)
xlabel('P2 value')
ylabel('Frequency')
xlim([0,1])
ylim([0,scalehist])

% pDif
Pdif_unlim = (P1_sim_unlim-P2_sim_unlim);
subplot(4,3,12)
histogram(Pdif_unlim,'BinWidth',.01)
xlabel('(P1 - P2) value')
ylabel('Frequency')
xlim([0,1])

%% Calculate the achieved power under the two scenarios above

p = .05; % What is our alpha?

% Where is the critical value between the two P1-P2 distributions?
crit = prctile(Pdif_unlim,(1-p)*100);
power_achieved = sum(Pdif_switch > crit)/length(Pdif_switch);
fprintf('Power acheived: %.1f%%\n',power_achieved*100)
ylim([0,scalehist])