% Chap3PlotMonoexpSignalWithKurtosisParameters.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap3PlotMonoexpSignalWithKurtosisParameters()

BValIncrement = 2;
BValueArray = 0:BValIncrement:30;
BDim = length(BValueArray);
STDNoise = 0.04;
RateDecay = 0.1;
KurtosisParam1 = -1;
KurtosisParam2 = 2;

NoiseFreeSignal = 1.*exp(-RateDecay.*BValueArray);
KurtosisSignal1 = 1.*exp(-RateDecay.*BValueArray + (KurtosisParam1.*RateDecay.^2.*BValueArray.^2)/6);
KurtosisSignal2 = 1.*exp(-RateDecay.*BValueArray + (KurtosisParam2.*RateDecay.^2.*BValueArray.^2)/6);
   
% Do a line plot of the true value    
figure('position',[50 50 800 600]);semilogy(BValueArray,NoiseFreeSignal,'Color','k','LineWidth',2,'LineStyle','--');hold all;ylim([0.01 1]);xlim([0 30]);
xlabel('b (AU)','fontsize',18)
ylabel('log Signal Amplitude','fontsize',18)
set(gca,'FontSize',14)

% Overlay the two
semilogy(BValueArray,KurtosisSignal1,'Color','b','LineWidth',2,'LineStyle','-');
semilogy(BValueArray,KurtosisSignal2,'Color','g','LineWidth',2,'LineStyle','-');

text(18,0.2,'Monoexponential Decay','Color','k','FontSize',16);
text(12,0.7,'Kurtosis = 2','Color','g','FontSize',16);
text(16,0.035,'Kurtosis = -1','Color','b','FontSize',16);