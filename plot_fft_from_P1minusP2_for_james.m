%% compute which frequencies will come out of the FFT analysis based on the
% sampling frequency and number of points
Fs = 1/.04; % sampling frequency
T = 1/Fs; % period
L = 13-1; % number of sampled periods
f = Fs*(0:(L/2))/L; freqs=f(2:end); % sampled frequencies

% here P1_all and P2_all should be 13*2*13 matrices with
%   dim1 = delays
%   dim2 = validity (1=invalid, 2=valid)
%   dim3 = observers
Pdiff=P1_all-P2_all;
fft_Pdiff=fft(Pdiff,13,1);
a_fft_Pdiff=abs(fft_Pdiff);

%%
plot_order=[2, 1];
figure('Position',get(groot,'ScreenSize'));
for val=1:2
    subplot(1,2,plot_order(val)); hold on;
    plot(f(2:end),mean(a_fft_Pdiff(2:size(f,2),val,:),3),'ko-',...
        'LineWidth',3,'MarkerFaceColor',[1 1 1],'MarkerSize',12,'Color',[.2 .2 .2])
    if val==2; title('Valid'); elseif val==1; title('Invalid'); end

    set(gca,'YTick',.5:.2:1.6,'FontSize',13,'LineWidth',2','Fontname','Ariel')
    set(gca,'XTick',f(2:end),'FontSize',13,'LineWidth',2','Fontname','Ariel')

    ylabel('Amplitude of oscillation in P1-P2','FontSize',12,'Fontname','Ariel')
    xlabel('Frequency (Hz)','FontSize',12,'Fontname','Ariel')
    ylim([.25 1.1]); xlim([f(2)-1 f(end)+1])
end
suptitle(sprintf('fft diff P1 minus P2 - %i subj - probeGratPos : %s ',n_obs, probeGratPos));
