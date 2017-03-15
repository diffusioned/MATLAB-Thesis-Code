function Chap2LinePlotOfSignalRange()

% This function is for creating a log-lin plot of the range of signals that will be tested
BValueArray = [0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6]; % from the PLOS ONE paper + 30
BDim = length(BValueArray);

% Maximum signal is D1 of 1 with a SF1 of 1 (Signal E)
MaxDecaySF = 1;
MaxDecayRate = 1;
MaxDecaySignal = MaxDecaySF.*exp(-MaxDecayRate.*BValueArray);

% Minimum decay signal is D2 of 1/20 (ratio of 20:1) with a SF1 of 0 (A2 = 1) (Signal A)
MinDecayA2 = 1;
MinDecayRate = 1/20.;
MinDecaySignal = MinDecayA2.*exp(-MinDecayRate.*BValueArray);


figure('Position', [50,50,1000,800]);
% Plot the noise floor values
semilogy(BValueArray,repmat(0.04,1,BDim),'LineStyle','--','LineWidth',2, 'Color','k');hold all; % SNR 25
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
xlabel('Diffusion Weighting (AU)','fontsize', 18, 'FontWeight', 'bold')
ylabel('log Signal Amplitude','fontsize', 18, 'FontWeight', 'bold')
% semilogy(BValueArray,repmat(0.01,1,BDim),'LineStyle','--','LineWidth',1.5, 'Color','c'); % SNR 100

% Add the thin b value lines
YRepValues = 0.01:0.005:1;
for i = 1:BDim
    semilogy([BValueArray(i), BValueArray(i)], [0.01, 1], 'LineWidth',1,'Color','r');
end

% Plot the min and max lines
semilogy(BValueArray,MaxDecaySignal,'LineWidth',2,'Color','b');ylim([0.01 1]);xlim([0 25.6]);  
semilogy(BValueArray,MinDecaySignal,'LineWidth',2,'Color','b');

% Add some text explanation for all added labels/points
text(15,0.6,'Minimum Decay Value','Color','b','FontSize',18, 'FontWeight', 'bold');
text(3,0.1,'Maximum Decay Value','Color','b','FontSize',18, 'FontWeight', 'bold');
text(18,0.035,'SNR 25 Noise Average','Color','k','FontSize',18 , 'FontWeight', 'bold');
