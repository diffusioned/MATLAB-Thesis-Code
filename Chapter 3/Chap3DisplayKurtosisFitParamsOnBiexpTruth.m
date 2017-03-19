% Chap3DisplayKurtosisFitParamsOnBiexpTruth.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap3DisplayKurtosisFitParamsOnBiexpTruth()

% Change the directories for displaying different data sets
% Load noise free parameter values for data organization
load('YourPath\NoiseFreeParameters.mat');

% Load kurtosis fits for SNR 25
load('YourPath\KurtFitArray.mat');

% Grab the dimensions
[NoiseFreeDim ~] = size(NFTestParameterArray);
[~, NoisySigDim, ~] = size(KurtFitArray);
TestParamDim = 3;  % Amp, ADC, K

% Sort the noise free parameters into bins - create 100x100 bins
SF1BinEdges = 0:0.01:1; SF1BinDim = 100;
D1D2BinEdges = 2:0.18:20; D1D2BinDim = 100;
EdgeVals = {SF1BinEdges, D1D2BinEdges}; % Cell array

NFSF1ValueArray = NFTestParameterArray(:,1)./(NFTestParameterArray(:,1)+NFTestParameterArray(:,2));
NFD1D2RatioArray = NFTestParameterArray(:,3)./NFTestParameterArray(:,4);
[SF1Cnt,~,SF1Idxs] = histcounts(NFSF1ValueArray,SF1BinEdges);
[D1D2Cnt,~,D1D2Idxs] = histcounts(NFD1D2RatioArray,D1D2BinEdges);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter Fits
% Loop through the bin indices and combine values by bin
BinStdDevArray = zeros(SF1BinDim, D1D2BinDim, TestParamDim);
BinCoeffVarArray = zeros(SF1BinDim, D1D2BinDim, TestParamDim);
BinIQRArray = zeros(SF1BinDim, D1D2BinDim, TestParamDim);
BinIQR_MedianArray = zeros(SF1BinDim, D1D2BinDim, TestParamDim);
AmpParameterEstimates =  squeeze(KurtFitArray(:,:,1));
ADCParameterEstimates =  squeeze(KurtFitArray(:,:,2));
KParameterEstimates =  squeeze(KurtFitArray(:,:,3));
% Combine all errors into distribution metrics by parameter
for i = 1:SF1BinDim
    for j = 1:D1D2BinDim
        CurBinIdxs = find(SF1Idxs == i & D1D2Idxs == j);
        AmpParamEst = AmpParameterEstimates(CurBinIdxs,:);
        ADCParamEst = ADCParameterEstimates(CurBinIdxs,:);
        KParamEst = KParameterEstimates(CurBinIdxs,:);

        % Standard deviation and dispersion measures compared to the fitted values not the deviations
        BinStdDevArray(i,j,1) = std(AmpParamEst(:));
        BinStdDevArray(i,j,2) = std(ADCParamEst(:));
        BinStdDevArray(i,j,3) = std(KParamEst(:));
        % CV
        BinCoeffVarArray(i,j,1) = std(AmpParamEst(:))./mean(AmpParamEst(:)); 
        BinCoeffVarArray(i,j,2) = std(ADCParamEst(:))./mean(ADCParamEst(:)); 
        BinCoeffVarArray(i,j,3) = std(KParamEst(:))./mean(KParamEst(:)); 
        % Robust dispersion method will be the interquartile range divided by the median
        BinIQRArray(i,j,1) = iqr(AmpParamEst(:));
        BinIQR_MedianArray(i,j,2) = iqr(ADCParamEst(:))./median(ADCParamEst(:)); 
        BinIQRArray(i,j,3) = iqr(KParamEst(:));
        
    end
end

% Fonts for plotting
PlotFontSize = 18;
ContourFontSize = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amp standard deviation
figure('Position', [0,0,800,800]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(BinStdDevArray(:,:,1))');set(h,'edgecolor','none');
% title('Std. Dev. \itAmp\rm\bf Estimates');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itADC_{1} \rm\bf/ \itADC_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% caxis([0.014 0.028]); % SNR 25
% caxis([0. 0.03]); % SNR 200
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'NorthOutside');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADC standard deviation
figure('Position', [0,0,800,800]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(BinStdDevArray(:,:,2))');set(h,'edgecolor','none');
% title('Std. Dev. \itADC\rm\bf Estimates');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itADC_{1} \rm\bf/ \itADC_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% caxis([0.0 0.09]); % SNR 25
% caxis([0. 0.15]); % SNR 200
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'NorthOutside');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K standard deviation
figure('Position', [0,0,1000,800]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(BinStdDevArray(:,:,3))');set(h,'edgecolor','none');
% title('Std. Dev. \itK\rm\bf Estimates');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% caxis([0. 0.7]); % SNR 25
% caxis([0. 0.15]); % SNR 200
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'EastOutside');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
