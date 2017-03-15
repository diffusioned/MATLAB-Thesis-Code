% Chap1ShowRicianBiasPlot.m
% MATLAB file for creating Figure 3 in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap1ShowRicianBiasPlot()

BValIncrement = 0.5;
BValueArray = 0:BValIncrement:5;
BDim = length(BValueArray);
STDNoise = 0.1;
RateDecay = 1;

NoisySignalArray = zeros(1,BDim);
NoiseFreeSignal = 1.*exp(-RateDecay.*BValueArray);

for b = 1:BDim
    NoisySignalArray(b) = AddMagnitudeNoise(NoiseFreeSignal(b), STDNoise);
end
   
% Do a line plot of the true value    
figure('position',[50 50 800 600]);plot(BValueArray,NoiseFreeSignal,'Color','k','LineWidth',1.5,'LineStyle',':');hold all;
xlabel('b-Value','fontsize',18)
ylabel('Signal Amplitude','fontsize',18)
set(gca,'FontSize',14)
scatter(BValueArray,NoisySignalArray,30,'ob','filled');

stophere = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = AddMagnitudeNoise(a, s)
m = sqrt((a + s*randn(1)).^2 + (s*randn(1)).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%