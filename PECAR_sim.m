% PECAR_sim.m
%
% First data set from Laura and Mehdi.
%
% Mehdi's instructions:
% It's a cell array with 13 entries (one for each subject) and in each 
% cell contains a [number_of_trials , 4] matrix where each line is a trial:
%
%   1st column: position of probe1 (from 1 to 10)
%   2nd column: position of probe2
%   3rd column: correct response for probe1 (boolean)
%   4th column: correct response for probe2
%
% Just as a reminder I've only sent you the valid ones.
% Also here is a figure of the possible positions and their "label" in the 
% matrix, in the figure gratings are presented but the probes could be 
% presented in any of the dotted circles (of course neither the numbers 
% nor the dotted circles were presented during the experiment).

clear all; close all

savestuff = 0;  % Do you want to save figures?

N = 1000; % Number of simulations (N = 1000, takes a long time)
nTrials = 1400; % Nunber of trials per simulation
scalehist = N/2;

datadr = 'C:\Users\Kit Moreland\Dropbox\UW\Research\DividedAttention\PECAR_DugueSenoussi\';
datafile = dir([datadr,'*.mat']);
load(fullfile(datadr, datafile.name))
data = probe_info_valid;

nsubj = length(data);

%% Replicate basic analysis
% Ignore probe position and just look at correct on pos 1 and 2

P1 = zeros(1,nsubj);
P2 = zeros(1,nsubj);

for ss = 1:nsubj
    nt = size(data{ss},1);  % Number of trials
    
    curd = data{ss}; % Pull current subj's data from cell array
    
    r1 = curd(:,3);
    r2 = curd(:,4);
    
    [P1(ss), P2(ss)] = quadsolve(r1, r2, nt);
    
    pc(ss,:) = mean([r1,r2]);
    
end

%%
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

if savestuff
    saveas(gcf,fullfile(datadr,'RealdataP1P2_valid.png'))
end

%% Basic simulation: Switching assumed
% Use the mean P1 and P2 and SEs to simulate the outcome and calculate
% power achieved.

m1 = mean(P1);
m2 = mean(P2);
SE1 = std(P1)/sqrt(nsubj);
SE2 = std(P2)/sqrt(nsubj);

P1_sim = nan(1,N);
P2_sim = nan(1,N);

for n = 1:N
    
    % Draw a value of P1 and P2 from the distribution estimated with the se
    p1_cur = SE1*randn(1) + m1;
    p2_cur = SE2*randn(1) + m2;
    
    resp = nan(nTrials,2);
    for i = 1:nTrials
        resp(i,1) = binornd(1,p1_cur);
        resp(i,2) = binornd(1,p2_cur);
        
        resp(i,:) = Shuffle(resp(i,:)); % Shuffle so that these responses could have been on either side
    end
    pc_sim(n,:) = mean(resp);

    [P1_sim(n), P2_sim(n)] = quadsolve(resp(:,1), resp(:,2), nTrials);
end

figure(1)
subplot(4,3,2)
hold on
plot(repmat([1,2],n,1)', [P1_sim;P2_sim], 'o-', 'Color', .5*[1,1,1])
errorbar([1,2], [mean(P1_sim),mean(P2_sim)], ...
    [std(P1_sim)/sqrt(N),std(P2_sim)/sqrt(N)], ...
    [std(P1_sim)/sqrt(N),std(P2_sim)/sqrt(N)], 'k','Linewidth',2)
xlim([.5,2.5])
ylim([0,1])
set(gca,'XTick',[1,2],'XTickLabel',{'P1','P2'},'Linewidth',1.5,'FontSize',14)
ylabel('Probability correct report')
title('Simulated data P1 and P2: Switching')

% P1
subplot(4,3,5)
histogram(P1_sim,'BinWidth',.05)
xlabel('P1 value')
ylabel('Frequency')
xlim([0,1])
ylim([0,scalehist])
% P2
subplot(4,3,8)
histogram(P2_sim,'BinWidth',.05)
xlabel('P2 value')
ylabel('Frequency')
xlim([0,1])
ylim([0,scalehist])

% pDif
Pdif_switch = (P1_sim-P2_sim);
subplot(4,3,11)
histogram(Pdif_switch,'BinWidth',.05)
xlabel('(P1 - P2) value')
ylabel('Frequency')
xlim([0,1])
ylim([0,scalehist])

%% Basic simulation: unlimited capacity assumed

% Now assume that performance on both sides is the same
m1 = mean(pc(:));
m2 = m1;
SE1 = std(pc(:))/sqrt(nsubj);
SE2 = SE1;

P1_sim2 = nan(1,N);
P2_sim2 = nan(1,N);

for n = 1:N
    
    % Draw a value of P1 and P2 from the distribution estimated with the se
    p1_cur = SE1*randn(1) + m1;
    p2_cur = SE2*randn(1) + m2;
    
    resp = nan(nTrials,2);
    for i = 1:nTrials
        resp(i,1) = binornd(1,p1_cur);
        resp(i,2) = binornd(1,p2_cur);
        
        resp(i,:) = Shuffle(resp(i,:)); % Shuffle so that these responses could have been on either side
    end
    pc_sim(n,:) = mean(resp);

    [P1_sim2(n), P2_sim2(n)] = quadsolve(resp(:,1), resp(:,2), nTrials);
end

figure(1)
subplot(4,3,3)
hold on
plot(repmat([1,2],n,1)', [P1_sim2;P2_sim2], 'o-', 'Color', .5*[1,1,1])
errorbar([1,2], [mean(P1_sim2),mean(P2_sim2)], ...
    [std(P1_sim2)/sqrt(N),std(P2_sim2)/sqrt(N)], ...
    [std(P1_sim2)/sqrt(N),std(P2_sim2)/sqrt(N)], 'k','Linewidth',2)
xlim([.5,2.5])
ylim([0,1])
set(gca,'XTick',[1,2],'XTickLabel',{'P1','P2'},'Linewidth',1.5,'FontSize',14)
ylabel('Probability correct report')
title('Simulated data P1 and P2: Unlimited')

% P1
subplot(4,3,6)
histogram(P1_sim2,'BinWidth',.05)
xlabel('P1 value')
ylabel('Frequency')
xlim([0,1])
ylim([0,scalehist])
% P2
subplot(4,3,9)
histogram(P2_sim2,'BinWidth',.05)
xlabel('P2 value')
ylabel('Frequency')
xlim([0,1])
ylim([0,scalehist])

% pDif
Pdif_unlim = (P1_sim2 - P2_sim2);
subplot(4,3,12)
histogram(Pdif_unlim,'BinWidth',.05)
xlabel('(P1 - P2) value')
ylabel('Frequency')
xlim([0,1])
ylim([0,scalehist])

%% Calculate the achieved power under the two scenarios above

p = .05; % What is our alpha?

% Where is the critical value between the two P1-P2 distributions?
crit = prctile(Pdif_unlim,(1-p)*100);
power_achieved = sum(Pdif_switch > crit)/length(Pdif_switch);
fprintf('Power acheived: %.1f%%\n',power_achieved*100)

