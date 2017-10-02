%% plot P1 and P2 computed from what I've sent james
data_loc='C:\Users\Kit Moreland\Dropbox\UW\Research\DividedAttention\PECAR_DugueSenoussi\'; % TO-CHANGE
val=2;
load([data_loc '13obs_probe_info_validity_and_delays_allsubjects_forjames'])
P1_allj=zeros(1,13); P2_allj=zeros(1,13);
pboth_allj=zeros(1,13); pnone_allj=zeros(1,13); pone_allj=zeros(1,13);
for obs_i=1:13
    pboth_allj(obs_i)=mean(sum(probe_info{obs_i}(validity_all{obs_i}==val,3:4),2)==2);
    pone_allj(obs_i)=mean(sum(probe_info{obs_i}(validity_all{obs_i}==val,3:4),2)==1); % actually not used
    pnone_allj(obs_i)=mean(sum(probe_info{obs_i}(validity_all{obs_i}==val,3:4),2)==0);
    [P1_allj(obs_i), P2_allj(obs_i)] = quadratic_analysis(pboth_allj(obs_i), pnone_allj(obs_i));
end

% plotting
figure; hold on;
for obs_i=1:13
    plot([1 2], [P1_allj(obs_i) P2_allj(obs_i)], 'ko-', 'Color',[.5 .5 .5]);
end
errorbar(1, mean(P1_allj), std(P1_allj)./sqrt(n_obs),...
    'ro-','LineWidth',3,'MarkerFaceColor',[1 1 1],'MarkerSize',12,'Color',[.96 .37 .15])
errorbar(2, mean(P2_allj), std(P2_allj)./sqrt(n_obs),...
    'go-','LineWidth',3,'MarkerFaceColor',[1 1 1],'MarkerSize',12,'Color',[.13 .7 .15])
xlim([.7, 2.5]); ylim([0,.81])