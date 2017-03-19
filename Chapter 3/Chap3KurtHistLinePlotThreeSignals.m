% Chap3KurtHistLinePlotThreeSignals.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap3KurtHistLinePlotThreeSignals()

% This function is for creating a log-lin plot of the range of signals that will be tested
BValueArray = [0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6]; % from the PLOS ONE paper + 30
BPlotValues = 0:0.05:25.6;
BDim = length(BValueArray);
ADC1 = 1.0;

% Biexp signal A (top center)
AmpA1 = 0.5;  AmpA2 = 0.5;
ADCA2 = ADC1/15.;
DecaySignalA = AmpA1.*exp(-ADC1.*BPlotValues) + AmpA2.*exp(-ADCA2.*BPlotValues);

% Biexp signal B (top right)
AmpB1 = 0.95;  AmpB2 = 0.05;
ADCB2 = ADC1/15.;
DecaySignalB = AmpB1.*exp(-ADC1.*BPlotValues) + AmpB2.*exp(-ADCB2.*BPlotValues);

% Biexp signal C (center)
AmpC1 = 0.05;  AmpC2 = 0.95;
ADCC2 = ADC1/15.;
DecaySignalC = AmpC1.*exp(-ADC1.*BPlotValues) + AmpC2.*exp(-ADCC2.*BPlotValues);

% Plot these signals
figure('Position', [50,0,1200,800]);
semilogy(BPlotValues,DecaySignalA,'LineWidth',4,'Color','b');ylim([0.01 1]);xlim([0 25.6]); hold on;
semilogy(BPlotValues,DecaySignalB,'LineWidth',4,'Color','g');
semilogy(BPlotValues,DecaySignalC,'LineWidth',4,'Color','r'); 

set(gca,'FontSize',16, 'FontWeight', 'bold');
xlabel('Diffusion Weighting (AU)', 'fontsize', 24, 'fontweight', 'bold');
ylabel('log Signal Amplitude', 'fontsize', 24, 'fontweight', 'bold');


% Add the thin b value lines
YRepValues = 0.01:0.005:1;
for i = 1:BDim
    semilogy([BValueArray(i), BValueArray(i)], [0.01, 1], 'LineWidth',1,'Color','r');
end

% Plot the SNR noise level
semilogy(BValueArray,repmat(0.04, 1, BDim),'LineWidth',2,'Color','k', 'LineStyle', '-.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add some text explanation for all added labels/points
text(16.5,0.035,'SNR 25 Noise Average','Color','k','FontSize',20 , 'FontWeight', 'bold');
text(14.5,0.14,'A','Color','b','FontSize',36, 'FontWeight', 'bold');
text(14.5,0.025,'B','Color','g','FontSize',36, 'FontWeight', 'bold');
text(14.5,0.45,'C','Color','r','FontSize',36, 'FontWeight', 'bold');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = AddMagnitudeNoise(a, s)
m = sqrt((a + s*randn(1)).^2 + (s*randn(1)).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%