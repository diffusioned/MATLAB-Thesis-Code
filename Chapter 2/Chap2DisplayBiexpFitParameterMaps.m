% Chap2DisplayBiexpFitParameterMaps.m
% MATLAB file for displaying fitted parameters in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap2DisplayBiexpFitParameterMaps()

% Change the directories for displaying different data sets
% Load noise free parameter values for data organization
load('YourPath\NoiseFreeParameters.mat');
% Load monoexp fits
load('YourPath\MonoexpFitArray.mat'); % SNR 25
% Load biexp fits
load('YourPath\BiexpFitArray.mat');

% TEST LOAD NOISY SIGNALS
% load('U:\Code\MATLAB\Thesis Chapter\Biexponential\Monte Carlo Data\Thesis 10M Noisy Signals\NoisyBiexpSignalsSNR200.mat');
% END TEST

% Get the dimensions
[NoiseFreeDim ParamDim] = size(NFTestParameterArray)
[~, NoisySigDim, ~] = size(BiexpFitArray);
TestParamDim = 3;  % For deviance, we just want SF1, D1, and D2
AbsParameterDevianceArray = zeros(NoiseFreeDim, NoisySigDim, TestParamDim);
PctParameterDevianceArray = zeros(NoiseFreeDim, NoisySigDim, TestParamDim);

% Loop through all noisy signals and calculate the deviance between measured parameter value and known true value for both absolute and percentage error
for s = 1:NoisySigDim
    % Calculate the deviance between the signal parameter and the known true value
    AbsParameterDevianceArray(:,s,1) = squeeze(BiexpFitArray(:,s,1))./(squeeze(BiexpFitArray(:,s,1))+squeeze(BiexpFitArray(:,s,2))) - ...
                                    squeeze(NFTestParameterArray(:,1))./(squeeze(NFTestParameterArray(:,1))+squeeze(NFTestParameterArray(:,2))); % SF1 abs error
    AbsParameterDevianceArray(:,s,2) = squeeze(BiexpFitArray(:,s,3)) - squeeze(NFTestParameterArray(:,3)); % ADC1 abs error
    AbsParameterDevianceArray(:,s,3) = squeeze(BiexpFitArray(:,s,4)) - squeeze(NFTestParameterArray(:,4)); % ADC2 abs error
    PctParameterDevianceArray(:,s,2) = (squeeze(BiexpFitArray(:,s,3)) - squeeze(NFTestParameterArray(:,3)))*100./squeeze(NFTestParameterArray(:,3)); % ADC1 pct error
    PctParameterDevianceArray(:,s,3) = (squeeze(BiexpFitArray(:,s,4)) - squeeze(NFTestParameterArray(:,4)))*100./squeeze(NFTestParameterArray(:,4)); % ADC2 pct error
end

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
BinMeanArray = zeros(SF1BinDim, D1D2BinDim, TestParamDim);
BinMedianArray = zeros(SF1BinDim, D1D2BinDim, TestParamDim);
BinStdDevArray = zeros(SF1BinDim, D1D2BinDim, TestParamDim);
BinCoeffVarArray = zeros(SF1BinDim, D1D2BinDim, TestParamDim);
BinIQRArray = zeros(SF1BinDim, D1D2BinDim, TestParamDim);
BinIQR_MedianArray = zeros(SF1BinDim, D1D2BinDim, TestParamDim);
SF1ParameterEstimates =  squeeze(BiexpFitArray(:,:,1))./(squeeze(BiexpFitArray(:,:,1))+squeeze(BiexpFitArray(:,:,2)));
ADC1ParameterEstimates =  squeeze(BiexpFitArray(:,:,3));
ADC2ParameterEstimates =  squeeze(BiexpFitArray(:,:,4));
BinMonoexpFitMeanArray = zeros(SF1BinDim, D1D2BinDim, 2); % Also monoexp fits
BinMonoexpFitStdDevArray = zeros(SF1BinDim, D1D2BinDim, 2); % Also monoexp fits
BinMonoADCCoeffVarArray = zeros(SF1BinDim, D1D2BinDim); % Also monoexp fits
% Combine all errors into distribution metrics by parameter
for i = 1:SF1BinDim
    for j = 1:D1D2BinDim
        CurBinIdxs = find(SF1Idxs == i & D1D2Idxs == j);
        AbsSF1DevArray = AbsParameterDevianceArray(CurBinIdxs,:,1);
        ADC1PctDevArray = PctParameterDevianceArray(CurBinIdxs,:,2);
        ADC2PctDevArray = PctParameterDevianceArray(CurBinIdxs,:,3);
        SF1ParamEst = SF1ParameterEstimates(CurBinIdxs,:);
        ADC1ParamEst = ADC1ParameterEstimates(CurBinIdxs,:);
        ADC2ParamEst = ADC2ParameterEstimates(CurBinIdxs,:);
        MonoAmpParamEst = MonoexpFitArray(CurBinIdxs,:,1);
        MonoADCParamEst = MonoexpFitArray(CurBinIdxs,:,2);

        BinMeanArray(i,j,1) = mean(AbsSF1DevArray(:)); %mean abs error
        BinMeanArray(i,j,2) = mean(ADC1PctDevArray(:)); % mean pct error
        BinMeanArray(i,j,3) = mean(ADC2PctDevArray(:)); % mean pct error
        BinMedianArray(i,j,1) = median(AbsSF1DevArray(:)); % median pct error
        BinMedianArray(i,j,2) = median(ADC1PctDevArray(:)); % median pct error
        BinMedianArray(i,j,3) = median(ADC2PctDevArray(:)); % median pct error
        % Standard deviation and dispersion measures compared to the fitted values not the deviations
        BinStdDevArray(i,j,1) = std(SF1ParamEst(:));
        BinCoeffVarArray(i,j,2) = std(ADC1ParamEst(:))./mean(ADC1ParamEst(:)); 
        BinCoeffVarArray(i,j,3) = std(ADC2ParamEst(:))./mean(ADC2ParamEst(:));
        % Robust dispersion method will be the interquartile range divided by the median
        BinIQRArray(i,j,1) = iqr(SF1ParamEst(:));
        BinIQR_MedianArray(i,j,2) = iqr(ADC1ParamEst(:))./median(ADC1ParamEst(:)); 
        BinIQR_MedianArray(i,j,3) = iqr(ADC2ParamEst(:))./median(ADC2ParamEst(:));
        BinMonoexpFitMeanArray(i,j,1) = mean(MonoAmpParamEst(:));
        BinMonoexpFitMeanArray(i,j,2) = mean(MonoADCParamEst(:));
        BinMonoexpFitStdDevArray(i,j,1) = std(MonoAmpParamEst(:));
        BinMonoexpFitStdDevArray(i,j,2) = std(MonoADCParamEst(:));
        BinMonoADCCoeffVarArray(i,j) = std(MonoADCParamEst(:))./mean(MonoADCParamEst(:));
        
    end
end

% Fonts for plotting
PlotFontSize = 18;
ContourFontSize = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SF1, mean absolute error and standard deviation
% figure('Position', [100,100,1000,800]); % For photoshop
figure('Position', [0,0,1000,800]); % For parameter value overlay
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(BinMeanArray(:,:,1))');set(h,'edgecolor','none');
title('Mean Absolute Error (color) and Std. Dev. (contours) of \itSF_{1}\rm\bf Estimates');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
caxis([-0.3 0.5]);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
colormap(jet); % gives a good contrast for negative numbers
% HandleC =  colorbar('location', 'SouthOutside');
% Overlay std dev contours
hold on;
Contours = 0.:0.05:0.4;
% The std dev contours need a bit of averaging for smoothness
flt_size = 3; avg_flt_size = 3; pad_size = flt_size;% window size
AvgFilter = fspecial('average', [avg_flt_size avg_flt_size]);
SF1StdDevPadded = padarray(BinStdDevArray(:,:,1)', [pad_size pad_size], 'replicate'); %zero gradient pad
SF1StdDev = imcrop(filter2(AvgFilter, SF1StdDevPadded), [pad_size + 1, pad_size + 1, SF1BinDim - 1, D1D2BinDim - 1]);
[Cons, HandleSF1Var] = contour(SF1BinEdges(2:end), D1D2BinEdges(2:end), SF1StdDev, Contours, 'Linecolor', [.0 .0 .0]);
set(HandleSF1Var, 'ShowText', 'on', 'TextStep', get(HandleSF1Var, 'LevelStep')*200., 'linewidth', 1.5, 'labelspacing', 300)
clabel(Cons,HandleSF1Var, 'FontSize', ContourFontSize);
pausehere = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADC1, mean pct error and coefficient of variation
figure('Position', [100,100,1000,800]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(BinMeanArray(:,:,2))');set(h,'edgecolor','none');
title('Mean Pct. Error (color) and Coeff. of Variation (contours) of \itADC_{1}\rm\bf Estimates');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itADC_{1} \rm\bf/ \itADC_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
caxis([-30 50]);
colormap(jet); % gives a good contrast for negative numbers
HandleC =  colorbar('location', 'SouthOutside');
% Overlay std dev contours
hold on;
Contours = 0.:0.1:1.2;
% The variance contours need a bit of filtering
flt_size = 3; avg_flt_size = 3; pad_size = flt_size;% window size
AvgFilter = fspecial('average', [avg_flt_size avg_flt_size]);
ADC1CoeffVarPadded = padarray(BinCoeffVarArray(:,:,2)', [pad_size pad_size], 'replicate'); %zero gradient pad
ADC1CoeffVar = imcrop(filter2(AvgFilter, ADC1CoeffVarPadded), [pad_size + 1, pad_size + 1, SF1BinDim - 1, D1D2BinDim - 1]);
[Cons, HandleADC1Var] = contour(SF1BinEdges(2:end), D1D2BinEdges(2:end), ADC1CoeffVar, Contours, 'Linecolor', [.0 .0 .0]);
set(HandleADC1Var, 'ShowText', 'on', 'TextStep', get(HandleADC1Var, 'LevelStep')*200., 'linewidth', 1.5, 'labelspacing', 300)
clabel(Cons,HandleADC1Var, 'FontSize', ContourFontSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADC2, mean pct error and coefficient of variation
figure('Position', [100,100,1000,800]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(BinMeanArray(:,:,3))');set(h,'edgecolor','none');
title('Mean Pct. Error (color) and Coeff. of Variation (contours) of \itADC_{2}\rm\bf Estimates');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itADC_{1} \rm\bf/ \itADC_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
caxis([-30 50]);
colormap(jet); % gives a good contrast for negative numbers
HandleC =  colorbar('location', 'SouthOutside');
% Overlay std dev contours
hold on;
Contours = 0.:0.1:1.2;
% The variance contours need a bit of filtering
flt_size = 3; avg_flt_size = 3; pad_size = flt_size;% window size
AvgFilter = fspecial('average', [avg_flt_size avg_flt_size]);
ADC2CoeffVarPadded = padarray(BinCoeffVarArray(:,:,3)', [pad_size pad_size], 'replicate'); %zero gradient pad
ADC2CoeffVar = imcrop(filter2(AvgFilter, ADC2CoeffVarPadded), [pad_size + 1, pad_size + 1, SF1BinDim - 1, D1D2BinDim - 1]);
[Cons, HandleADC2Var] = contour(SF1BinEdges(2:end), D1D2BinEdges(2:end), ADC2CoeffVar, Contours, 'Linecolor', [.0 .0 .0]);
set(HandleADC2Var, 'ShowText', 'on', 'TextStep', get(HandleADC2Var, 'LevelStep')*200., 'linewidth', 1.5, 'labelspacing', 300)
clabel(Cons,HandleADC2Var, 'FontSize', ContourFontSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REPEAT FOR ROBUST STATISTICS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SF1, median absolute error and IQR
figure('Position', [100,100,1000,800]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(BinMedianArray(:,:,1))');set(h,'edgecolor','none');
title('Median Absolute Error (color) and IQR (contours) of \itSF_{1}\rm\bf Estimates');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itADC_{1} \rm\bf/ \itADC_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
caxis([-0.3 0.5]);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
colormap(jet); % gives a good contrast for negative numbers
HandleC =  colorbar('location', 'SouthOutside');
% Overlay std dev contours
hold on;
Contours = 0.1:0.1:0.8;
% The std dev contours need a bit of averaging for smoothness
flt_size = 3; avg_flt_size = 3; pad_size = flt_size;% window size
AvgFilter = fspecial('average', [avg_flt_size avg_flt_size]);
SF1IQRPadded = padarray(BinIQRArray(:,:,1)', [pad_size pad_size], 'replicate'); %zero gradient pad
SF1IQR = imcrop(filter2(AvgFilter, SF1IQRPadded), [pad_size + 1, pad_size + 1, SF1BinDim - 1, D1D2BinDim - 1]);
[Cons, HandleSF1Var] = contour(SF1BinEdges(2:end), D1D2BinEdges(2:end), SF1IQR, Contours, 'Linecolor', [.0 .0 .0]);
set(HandleSF1Var, 'ShowText', 'on', 'TextStep', get(HandleSF1Var, 'LevelStep')*200., 'linewidth', 1, 'labelspacing', 300)
clabel(Cons,HandleSF1Var, 'FontSize', ContourFontSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADC1, median pct error and IQR/median
figure('Position', [100,100,1000,800]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(BinMedianArray(:,:,2))');set(h,'edgecolor','none');
title('Median Pct. Error (color) and IQR/Median (contours) of \itADC_{1}\rm\bf Estimates');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itADC_{1} \rm\bf/ \itADC_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
caxis([-60 20]);
colormap(jet); % gives a good contrast for negative numbers
HandleC =  colorbar('location', 'SouthOutside');
% Overlay std dev contours
hold on;
Contours = 0.:0.1:1.5;
% The variance contours need a bit of filtering
flt_size = 3; avg_flt_size = 3; pad_size = flt_size;% window size
AvgFilter = fspecial('average', [avg_flt_size avg_flt_size]);
ADC1IQRMedianPadded = padarray(BinIQR_MedianArray(:,:,2)', [pad_size pad_size], 'replicate'); %zero gradient pad
ADC1IQRMedian = imcrop(filter2(AvgFilter, ADC1IQRMedianPadded), [pad_size + 1, pad_size + 1, SF1BinDim - 1, D1D2BinDim - 1]);
[Cons, HandleADC1IQRMedian] = contour(SF1BinEdges(2:end), D1D2BinEdges(2:end), ADC1IQRMedian, Contours, 'Linecolor', [.0 .0 .0]);
set(HandleADC1IQRMedian, 'ShowText', 'on', 'TextStep', get(HandleADC1IQRMedian, 'LevelStep')*200., 'linewidth', 1, 'labelspacing', 300)
clabel(Cons,HandleADC1IQRMedian, 'FontSize', ContourFontSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADC2, median pct error and IQR/median
figure('Position', [100,100,1000,800]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(BinMedianArray(:,:,3))');set(h,'edgecolor','none');
title('Median Pct. Error (color) and IQR/Median (contours) of \itADC_{2}\rm\bf Estimates');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itADC_{1} \rm\bf/ \itADC_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
caxis([-60 20]);
colormap(jet); % gives a good contrast for negative numbers
HandleC =  colorbar('location', 'SouthOutside');
% Overlay std dev contours
hold on;
Contours = 0.:0.1:1.5;
% The variance contours need a bit of filtering
flt_size = 3; avg_flt_size = 3; pad_size = flt_size;% window size
AvgFilter = fspecial('average', [avg_flt_size avg_flt_size]);
ADC2IQRMedianPadded = padarray(BinIQR_MedianArray(:,:,3)', [pad_size pad_size], 'replicate'); %zero gradient pad
ADC2IQRMedian = imcrop(filter2(AvgFilter, ADC2IQRMedianPadded), [pad_size + 1, pad_size + 1, SF1BinDim - 1, D1D2BinDim - 1]);
[Cons, HandleADC2IQRMedian] = contour(SF1BinEdges(2:end), D1D2BinEdges(2:end), ADC2IQRMedian, Contours, 'Linecolor', [.0 .0 .0]);
set(HandleADC2IQRMedian, 'ShowText', 'on', 'TextStep', get(HandleADC2IQRMedian, 'LevelStep')*200., 'linewidth', 1, 'labelspacing', 300)
clabel(Cons,HandleADC2IQRMedian, 'FontSize', ContourFontSize);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD IN MONOEXPONENTIAL PARAMETER ESTIMATES
% No true value for the monoexponential, so compare the variation of the monoexponential estimates at all points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amplitude colormap for standard deviation
fS0 = figure('Position', [0,0,800,600]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(BinMonoexpFitStdDevArray(:,:,1))');set(h,'edgecolor','none');
title('Standard Deviation of \itS_{0}\rm\bf Estimates');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
caxis([0.014 0.028]);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
% HandleC =  colorbar('location', 'NorthOutside');

% ADC colormap for standard deviation
figure('Position', [100,100,1000,800]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(BinMonoexpFitStdDevArray(:,:,2))');set(h,'edgecolor','none');
title('Standard Deviation of \itADC\rm\bf Estimates');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itADC_{1} \rm\bf/ \itADC_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'NorthOutside');

% ADC colormap for mean
figure('Position', [100,100,1000,800]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(BinMonoexpFitMeanArray(:,:,2))');set(h,'edgecolor','none');
title('Mean of \itADC\rm\bf Estimates');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itADC_{1} \rm\bf/ \itADC_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'SouthOutside');

% ADC colormap for coefficient of variation
figure('Position', [0,0,800,600]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(BinMonoADCCoeffVarArray(:,:))');set(h,'edgecolor','none');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
caxis([0.07 0.16]);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC = colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'CV');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];

