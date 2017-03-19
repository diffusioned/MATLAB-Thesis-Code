% Chap3DisplayKurtosisFitParamsOnMonoexpTruth.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap3DisplayKurtosisFitParamsOnMonoexpTruth()

% Load noise free parameter values for data organization
load('YourPath\MonoexpNFParametersForKurtosis.mat');

% Load kurtosis RSS values
load('YourPath\KurtResidualArray11BValuesKurtChap.mat');

% Load kurtosis fits
load('YourPath\KurtFitArray11BValuesKurtChap.mat');

% Load the Jacobian condition
load('YourPath\KurtJacCondToMonoArray.mat');


% TEST LOAD NOISY SIGNALS
% load('U:\Code\MATLAB\Thesis Chapter\Biexponential\Monte Carlo Data\Thesis 10M Noisy Signals\NoisyBiexpSignalsSNR200.mat');
% END TEST

% Grab the dimensions
[NoiseFreeDim ~] = size(NFTestParameterArray);
[~, NoisySigDim, ~] = size(KurtFitArray);
TestParamDim = 3;  % Amp, ADC, K

% Sort the noise free parameters into bins - create 96 bins
ADCBinEdges = 0.05:0.01:1.; ADCBinDim = 95;
EdgeVals = {ADCBinEdges}; % Cell array
[ADCCnt,~,ADCIdxs] = histcounts(NFTestParameterArray(:,2),ADCBinEdges);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Fits
% Loop through the bin indices and combine values by bin
BinStdDevArray = zeros(ADCBinDim, TestParamDim);
BinCoeffVarArray = zeros(ADCBinDim, TestParamDim);
BinIQRArray = zeros(ADCBinDim, TestParamDim);
BinIQR_MedianArray = zeros(ADCBinDim, TestParamDim);
AmpParameterEstimates =  squeeze(KurtFitArray(:,:,1));
ADCParameterEstimates =  squeeze(KurtFitArray(:,:,2));
KParameterEstimates =  squeeze(KurtFitArray(:,:,3));
BinMeanJacCondArray = zeros(ADCBinDim,1);
% Combine all errors into distribution metrics by parameter
for i = 1:ADCBinDim
    CurBinIdxs = find(ADCIdxs == i);
    AmpParamEst = AmpParameterEstimates(CurBinIdxs,:);
    ADCParamEst = ADCParameterEstimates(CurBinIdxs,:);
    KParamEst = KParameterEstimates(CurBinIdxs,:);
    JacCondEst = KurtJacCondArray(CurBinIdxs,:);

    % Standard deviation and dispersion measures compared to the fitted values not the deviations
    BinStdDevArray(i,1) = std(AmpParamEst(:));
    BinStdDevArray(i,2) = std(ADCParamEst(:));
    BinStdDevArray(i,3) = std(KParamEst(:));
    % CV
    BinCoeffVarArray(i,1) = std(AmpParamEst(:))./mean(AmpParamEst(:)); 
    BinCoeffVarArray(i,2) = std(ADCParamEst(:))./mean(ADCParamEst(:)); 
    BinCoeffVarArray(i,3) = std(KParamEst(:))./mean(KParamEst(:)); 
    % Robust dispersion method will be the interquartile range divided by the median
    BinIQRArray(i,1) = iqr(AmpParamEst(:));
    BinIQR_MedianArray(i,2) = iqr(ADCParamEst(:))./median(ADCParamEst(:)); 
    BinIQRArray(i,3) = iqr(KParamEst(:));
    % Condition
    BinMeanJacCondArray(i) = mean(JacCondEst(:));
end

% Fonts for plotting
PlotFontSize = 14;
ContourFontSize = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amp standard deviation
figure('Position', [0,0,1000,800]);
subplot(2,2,1);
h1=bar(ADCBinEdges(2:end),squeeze(BinStdDevArray(:,1))','histc');
% title('Std. Dev. \itAmp\rm\bf Estimates');
xlabel('True \itADC\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Std. Dev. \itS\rm\bf(0) Estimates','FontWeight','bold','fontsize', PlotFontSize);
xlim([0.06 1.0]);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
% HandleC =  colorbar('location', 'NorthOutside');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADC standard deviation
% figure('Position', [0,0,1000,800]);
subplot(2,2,2);
h2=bar(ADCBinEdges(2:end),squeeze(BinStdDevArray(:,2))','histc');
% title('Std. Dev. \itADC\rm\bf Estimates');
xlabel('True \itADC\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Std. Dev. \itD_{App}\rm\bf Estimates','FontWeight','bold','fontsize', PlotFontSize);
xlim([0.06 1.0]);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K standard deviation
% figure('Position', [0,0,1000,800]);
subplot(2,2,3);
h3=bar(ADCBinEdges(2:end),squeeze(BinStdDevArray(:,3))','histc');
% title('Std. Dev. \itK\rm\bf Estimates');
xlabel('True \itADC\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Std. Dev. \itK_{App}\rm\bf Estimates','FontWeight','bold','fontsize', PlotFontSize);
xlim([0.06 1.0]);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean Jacobian Condition
subplot(2,2,4);
h3=bar(ADCBinEdges(2:end),squeeze(BinMeanJacCondArray(:))','histc');
% title('Std. Dev. \itK\rm\bf Estimates');
xlabel('True \itADC\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Mean Jacobian Condition','FontWeight','bold','fontsize', PlotFontSize);
xlim([0.06 1.0]);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%