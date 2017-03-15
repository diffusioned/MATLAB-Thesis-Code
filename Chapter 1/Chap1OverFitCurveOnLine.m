% Chap1ShowRicianBiasPlot.m
% MATLAB file for creating Figure 5 in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap1OverFitCurveOnLine()

BValIncrement = 0.5;
BValueArray = 0:BValIncrement:5;
BDim = length(BValueArray);
STDNoise = 0.1;
RateDecay = 1;

NoisySignalArray = zeros(1,BDim);
NoiseFreeSignal = 1.*exp(-RateDecay.*BValueArray);

for b = 1:BDim
    NoisySignalArray(b) = AddGaussianNoise(NoiseFreeSignal(b), STDNoise);
end
   
% Do a line plot of the true value    
figure('position',[50 50 800 600]);plot(BValueArray,NoiseFreeSignal,'Color','k','LineWidth',1.5,'LineStyle',':');hold all;
xlabel('b-Value','fontsize',18)
ylabel('Signal Amplitude','fontsize',18)
set(gca,'FontSize',14)
scatter(BValueArray,NoisySignalArray,30,'ob','filled');

% Fit the data
pRet1 = polyfit(BValueArray,NoisySignalArray, 1); % 1 means linear fit
pRet3 = polyfit(BValueArray,NoisySignalArray, 3); % 3 means cubic fit
LinearFit = pRet1(1)*BValueArray + pRet1(2);
Poly3Fit = pRet3(1)*(BValueArray.^3) + pRet3(2)*(BValueArray.^2) + pRet3(3)*BValueArray + pRet3(4);

plot(BValueArray,LinearFit,'Color','r','LineWidth',1.5,'LineStyle','-');
plot(BValueArray,Poly3Fit,'Color','g','LineWidth',1.5,'LineStyle','-');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = AddGaussianNoise(a, s)
m = a + s*randn(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
