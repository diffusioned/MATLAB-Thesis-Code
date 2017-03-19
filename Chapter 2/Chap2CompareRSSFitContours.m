% Chap2CompareRSSFitContours.m
% MATLAB file for displaying contours in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap2CompareRSSFitContours()
% Calculate the noise std. dev.
STDNoise = 0.04 % For SNR of 25 = 0.04

% Test Parameters
SF1 = 0.5;
ADC1 = 1.0; % Normalized to one
D1D2Ratio = 20.; % Choose 1, 2, or 20 for the plots below.
ADC2 = ADC1/D1D2Ratio;
Amp0 = 1.0;

BValueArray = [0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6]; % from the PLOS ONE paper
[BDim BDim] = size(BValueArray);
SignalDim = BDim; % Number of values to use in this test (or set to max = BDim)

% Add noise to the noise free signal
NoiseFreeSignal = CreateNoiseFreeSignal(Amp0, SF1, ADC1, ADC2, BValueArray);

% GaussianNoiseSignal = zeros(1,BDim);
RicianNoiseSignal = zeros(1,BDim);

% Add Noise to signal
for b = 1:length(BValueArray)
%     GaussianNoiseSignal(b) = AddGaussianNoise(NoiseFreeSignal(b), STDNoise);
    RicianNoiseSignal(b) = AddRicianNoise(NoiseFreeSignal(b), STDNoise);
end

% Create parameter grid to test the map over
SF1ParameterRegion = 0:0.005:1; SF1Dim = length(SF1ParameterRegion);
ADC1ParameterRegion = 0:0.02:4; ADC1Dim = length(ADC1ParameterRegion);
% Uncomment the value below you wish to display
ADC2ParameterRegion = 0:0.0025:0.5; ADC2Dim = length(ADC2ParameterRegion); % D1/D2 = 2 or 20
% ADC2ParameterRegion = 0.9:0.01:1.1; ADC2Dim = length(ADC2ParameterRegion); % D1/D2 = 1
% Save the RSS values in this array
RSSToNoiseFreeCalculationArray = zeros(SF1Dim,ADC1Dim,ADC2Dim);

% Now compare the parametric signal at each point with the initial signal and save the RSS
for i = 1:SF1Dim
    for j = 1:ADC1Dim
        for k = 1:ADC2Dim
            % Create the signal at this value
            CurrentSignal = CreateNoiseFreeSignal(Amp0, SF1ParameterRegion(i), ADC1ParameterRegion(j), ADC2ParameterRegion(k), BValueArray);
            % Loop through and find the sum of squares of the residuals of this signal compared to the noisy signal.
            RSSToNoiseFreeCalculationArray(i,j,k) = CompareRSSOfTwoSignals(CurrentSignal, NoiseFreeSignal);
        end
    end
    i
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PlotFontSize = 18;
% For 1/20 signal
figure('Position', [0,0,800,600]);
contour(SF1ParameterRegion, ADC1ParameterRegion, squeeze(RSSToNoiseFreeCalculationArray(:,:,21)'), 0:0.0005:0.1);
title('RSS Contour Map for \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio = 20');
xlabel('\itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('\itD_{1} \rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'RSS');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % For 1/2 signal
% figure('Position', [0,0,800,600]);
% contour(SF1ParameterRegion, ADC1ParameterRegion, squeeze(RSSToNoiseFreeCalculationArray(:,:,201)'), 0:0.0005:0.1);
% title('RSS Contour Map for \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio = 2');
% xlabel('\itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% ylabel('\itD_{1} \rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
% HandleC =  colorbar('location', 'EastOutside');
% HandleCLabel = ylabel(HandleC,'RSS');
% drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
% LabelPos = HandleCLabel.Position;
% HandleCLabel.Rotation = -90;
% HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % For 1/1 signal
% figure('Position', [0,0,800,600]);
% contour(SF1ParameterRegion, ADC1ParameterRegion, squeeze(RSSToNoiseFreeCalculationArray(:,:,11)'), 0:0.0005:0.1);
% title('RSS Contour Map for \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio = 1');
% xlabel('\itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% ylabel('\itD_{1} \rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
% HandleC =  colorbar('location', 'EastOutside');
% HandleCLabel = ylabel(HandleC,'RSS');
% drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
% LabelPos = HandleCLabel.Position;
% HandleCLabel.Rotation = -90;
% HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the minimum value in each
[minval, idx] = min(RSSToNoiseFreeCalculationArray(:))
[SF1Min, ADC1Min, ADC2Min] = ind2sub(size(RSSToNoiseFreeCalculationArray), idx)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = CreateNoiseFreeSignal(a0, sf1, adc1, adc2, b)
a = a0.*(sf1.*exp(-adc1.*b) + (1. - sf1).*exp(-adc2.*b));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = AddRicianNoise(a, s)
m = sqrt((a + s*randn(1))^2 + (s*randn(1))^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RSSReturn = CompareRSSOfTwoSignals(Signal1, Signal2)
RSSReturn = sum((Signal1 - Signal2).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%