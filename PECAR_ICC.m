% PECAR_ICC.m
%
% This is an intraclass cluster correlation (ICC) analysis to look for subject differences. If there are
% differences then we might want to include these in our models.
%
% 9/25/2017 JC Moreland

clear all; close all

datadr = 'C:\Users\Kit Moreland\Dropbox\UW\Research\DividedAttention\PECAR_DugueSenoussi\';
datafile = dir([datadr,'*.mat']);
load(fullfile(datadr, datafile.name))
data = probe_info_valid;

%% Open each observer and get a within subject variance
nsubj = length(data);

for ss = 1:nsubj
    nt = size(data{ss},1);  % Number of trials
    
    curd = data{ss}; % Pull current subj's data from cell array
    
    r1 = curd(:,3);
    r2 = curd(:,4);
    
    subjmeans(ss) = mean([r1;r2]);
    
    wVss(ss) = var([r1;r2]);
    
end

bV = var(subjmeans);        % 
wV = mean(wVss);    % Taking a mean of the subject variances

ICC = bV./(bV+wV);
fprintf(1,'ICC = %f\n',ICC)