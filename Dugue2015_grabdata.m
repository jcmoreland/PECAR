% Dugue2015_grabdata.m
%
% Use the data from figure 3 of Dugue et al 2015 to explore power. Note
% this is not strictly power acheived because I don't know the summary
% statistics of the inidividual subjects.
%
% Data was taken from a figure that used 15 SOAs. I use these 15 data
% points just to help estimate the standard error. A fun
% alternative would be to use each point like a new observer or to actually
% recreate the figure more.

clear all
close all

datadr = pwd;

%% Load data from grabit

load(fullfile(datadr,'Dugue2015_fig3_p1.mat'))
load(fullfile(datadr,'Dugue2015_fig3_p1se.mat'))
load(fullfile(datadr,'Dugue2015_fig3_p2.mat'))
load(fullfile(datadr,'Dugue2015_fig3_p2se.mat'))

p1_Dug = Dugue2015_fig3_p1(:,2);
p1_Dugse = Dugue2015_fig3_p1se(:,2) - p1_Dug;
p2_Dug = Dugue2015_fig3_p2(:,2);
p2_Dugse = Dugue2015_fig3_p2se(:,2) - p2_Dug;

t = linspace(30, 450, 15); % These are the actual x axis coordinates

%% recreate the Dugue et al 2015 paper figure
figure(1)
clf
hold on

errorbar(t, p1_Dug, p1_Dugse, p1_Dugse,'LineWidth',1.5)
errorbar(t, p2_Dug, p2_Dugse, p2_Dugse,'LineWidth',1.5)
set(gca,'LineWidth',1.5)
ylabel('Probe report probability')
xlabel('Time from search array onset (ms)')
title('P1 versus P2')

%% Use these p1/p2 standard errors for the null hypothesis that attention can be equally divided

p1_uc = .6; % unlimited capacity model for divided attention where performance is the same on both sides
p2_uc = p1_uc;
p1_se = mean(p1_Dugse);
p2_se = mean(p2_Dugse);

% I think the closed form solution would still be fine here but I will
% simulate just in case.

% nTrials = 96;   % Dugue 2015 report 1440 trials per subject but I think this was over the 15 SOAs fo this is 96 trials per sucject.
nTrials = 500;
N = 1000;       % Number of simulations

P1 = nan(N,1);
P2 = nan(N,1);

for n = 1:N
    
    % Draw a value of P1 and P2 from the distribution estimated with the se
    p1_cur = p1_se*randn(1) + p1_uc;
    p2_cur = p2_se*randn(1) + p2_uc;
    
    resp = nan(nTrials,2);
    for i = 1:nTrials
        resp(i,1) = binornd(1,p1_cur);
        resp(i,2) = binornd(1,p2_cur);
        
        resp(i,:) = Shuffle(resp(i,:)); % Shuffle so that these responses could have been on either side
        
    end

    Pboth = sum(resp(:,1) == 1 & resp(:,2)== 1)/nTrials;
    Pnone = sum(resp(:,1) == 0 & resp(:,2)== 0)/nTrials;
    
    b = 1 + Pboth - Pnone;
    c = Pboth;
    
    det = b^2 - 4*c;
    
    tmp = [(b + sign(det)*sqrt(abs(det)))/2,(b - sign(det)*sqrt(abs(det)))/2];
    
    P1(n) = max(tmp);
    P2(n) = min(tmp);
    
end

pc = mean(resp);
Pdif_unlim = abs(P1-P2);  % Take the absolute because the larger value on any analysis was taken as the most attended location

%% Figure for the simulation results 
figure(3)
clf
hold on

subplot(4,2,1)
errorbar([1,2],pc,[std(pc),std(pc)],[std(pc),std(pc)],'o')
xlim([0,3])
ylim([.25,1])
ylabel('Proportion correct')
title('Unlimited Capacity Performance')
set(gca,'XTick',[1,2],'XTickLabel',{'Left','Right'})

subplot(4,2,3)
hold on
histogram(P1,'BinWidth',0.05)
title('P1 probability correct report')
xlim([0,1])
subplot(4,2,5)
histogram(P2,'BinWidth',0.05)
title('P2 probability correct report')
xlim([0,1])

subplot(4,2,7)
histogram(Pdif_unlim,'BinWidth',0.05)
ylabel('Number of simulations')
xlabel('P1-P2')
xlim([-1,1])

%% Use these p1/p2 standard errors for the alternative hypothesis that attention oscillates

p1_uc = .6; % all or none switching model says that if you are right on one side you guess on the other
p2_uc = 1/6;
p1_se = mean(p1_Dugse);
p2_se = mean(p2_Dugse);

P1 = nan(N,1);
P2 = nan(N,1);

for n = 1:N
    
    % Draw a value of P1 and P2 from the distribution estimated with the se
    p1_cur = p1_se*randn(1) + p1_uc;
    p2_cur = p2_se*randn(1) + p2_uc;
    
    resp = nan(nTrials,2);
    for i = 1:nTrials
        resp(i,1) = binornd(1,p1_cur);
        resp(i,2) = binornd(1,p2_cur);
        
        resp(i,:) = Shuffle(resp(i,:)); % Shuffle so that these responses could have been on either side
    end
    
    Pboth = sum(resp(:,1) == 1 & resp(:,2)== 1)/nTrials;
    Pnone = sum(resp(:,1) == 0 & resp(:,2)== 0)/nTrials;
    
    b = 1 + Pboth - Pnone;
    c = Pboth;
    
    det = b^2 - 4*c;
    
    tmp = [(b + sign(det)*sqrt(abs(det)))/2,(b - sign(det)*sqrt(abs(det)))/2];
    
    P1(n) = max(tmp);
    P2(n) = min(tmp);
end

pc = mean(resp); % Shuffle because left and right are kind of meaningless here
Pdif_switch = abs(P1-P2);  % Take the absolute because the larger value on any analysis was taken as the most attended location

%% Figure for the simulation results 

subplot(4,2,2)
errorbar([1,2],pc,[std(pc),std(pc)],[std(pc),std(pc)],'o')
xlim([0,3])
ylim([.25,1])
ylabel('Proportion correct')
title('Switching Performance')
set(gca,'XTick',[1,2],'XTickLabel',{'Left','Right'})

subplot(4,2,4)
hold on
histogram(P1,'BinWidth',0.05)
xlim([0,1])
title('P1 probability correct report')
subplot(4,2,6)
histogram(P2,'BinWidth',0.05)
title('P2 probability correct report')
xlim([0,1])

subplot(4,2,8)
histogram(Pdif_switch,'BinWidth',0.05)
ylabel('Number of simulations')
xlabel('P1-P2')
xlim([-1,1])

%% Calculate the achieved power under the two scenarios above

p = .05; % What is our alpha

% Where is the critical value between the two P1-P2 distributions?
crit = prctile(Pdif_unlim,(1-p)*100);
power_achieved = sum(Pdif_switch > crit)/length(Pdif_switch)
