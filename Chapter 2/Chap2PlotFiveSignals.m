% Chap2PlotFiveSignals.m
% MATLAB file for plotting five signals in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap2PlotFiveSignals()

% This function is for creating a log-lin plot of the range of signals that will be tested
BValueArray = [0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6]; % from the PLOS ONE paper + 30
BPlotValues = 0:0.05:25.6;
BDim = length(BValueArray);
ADC1 = 1.0;

% Biexp signal A is near minimum value (top left)
AmpA1 = 0.05;  AmpA2 = 0.95;
ADCA2 = ADC1/19.;
DecaySignalA = AmpA1.*exp(-ADC1.*BPlotValues) + AmpA2.*exp(-ADCA2.*BPlotValues);

% Biexp signal B (top right)
AmpB1 = 0.95;  AmpB2 = 0.05;
ADCB2 = ADC1/19.;
DecaySignalB = AmpB1.*exp(-ADC1.*BPlotValues) + AmpB2.*exp(-ADCB2.*BPlotValues);

% Biexp signal C (center)
AmpC1 = 0.5;  AmpC2 = 0.5;
ADCC2 = ADC1/11.;
DecaySignalC = AmpC1.*exp(-ADC1.*BPlotValues) + AmpC2.*exp(-ADCC2.*BPlotValues);

% Biexp signal D (bottom left)
AmpD1 = 0.05;  AmpD2 = 0.95;
ADCD2 = ADC1/3.;
DecaySignalD = AmpD1.*exp(-ADC1.*BPlotValues) + AmpD2.*exp(-ADCD2.*BPlotValues);

% Biexp signal E is near maximum value (bottom right)
AmpE1 = 0.95;  AmpE2 = 0.05;
ADCE2 = ADC1/3.;
DecaySignalE = AmpE1.*exp(-ADC1.*BPlotValues) + AmpE2.*exp(-ADCE2.*BPlotValues);

% Plot these signals
figure('Position', [50,0,1000,1000]);
semilogy(BPlotValues,DecaySignalA,'LineWidth',4,'Color','b');ylim([0.01 1]);xlim([0 25.6]); hold on;
semilogy(BPlotValues,DecaySignalB,'LineWidth',4,'Color','g');
semilogy(BPlotValues,DecaySignalC,'LineWidth',4,'Color','r'); 
semilogy(BPlotValues,DecaySignalD,'LineWidth',4,'Color','c');
semilogy(BPlotValues,DecaySignalE,'LineWidth',4,'Color','m');

set(gca,'FontSize',16, 'FontWeight', 'bold');
xlabel('Diffusion Weighting (AU)', 'fontsize', 20, 'fontweight', 'bold');
ylabel('log Signal Amplitude', 'fontsize', 20, 'fontweight', 'bold');

% Add the thin b value lines
YRepValues = 0.01:0.005:1;
for i = 1:BDim
    semilogy([BValueArray(i), BValueArray(i)], [0.01, 1], 'LineWidth',1,'Color','r');
end

% Plot the SNR noise level
semilogy(BValueArray,repmat(0.04, 1, BDim),'LineWidth',2,'Color','k', 'LineStyle', '-.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add some text explanation for all added labels/points
% text(15,0.6,'Minimum Decay Value','Color','b','FontSize',18);
% text(3,0.1,'Maximum Decay Value','Color','b','FontSize',18);
text(18,0.035,'SNR 25 Noise Average','Color','k','FontSize',20 , 'FontWeight', 'bold');
text(15,0.55,'A','Color','b','FontSize',24, 'FontWeight', 'bold');
text(15,0.0288,'B','Color','g','FontSize',24, 'FontWeight', 'bold');
text(15,0.16,'C','Color','r','FontSize',24, 'FontWeight', 'bold');
text(8.2,0.08,'D','Color','c','FontSize',24, 'FontWeight', 'bold');
text(5,0.022,'E','Color','m','FontSize',24, 'FontWeight', 'bold');
