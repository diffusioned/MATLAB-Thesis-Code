% Chap2DisplayBootstrapBiexpFitMaps.m
% MATLAB file for display bootstrap data histograms in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap2DisplayBootstrapBiexpFitMaps()

% Change the directories for displaying different data sets
% Load noise free parameter values for data organization
load('YourPath\BootstrapNoiseFreeParameters.mat');

% Load biexp fits
load('YourPath\BiexpFitArray.mat');
load('YourPath\BiexpBootFitArray.mat'); % The 1000 bootstrap samples for each noisy signal
load('YourPath\BiexpFitJacobianArray.mat');
load('YourPath\BiexpResidualArray.mat');

% Get the dimensions
[NoiseFreeDim ParamDim] = size(BootstrapNFTestParameterArray)
[~, NoisySigDim, ~] = size(BiexpFitArray);
TestParamDim = 3;  % For deviance, we just want SF1, D1, and D2
AbsParameterDevianceArray = zeros(NoiseFreeDim, NoisySigDim, TestParamDim);
PctParameterDevianceArray = zeros(NoiseFreeDim, NoisySigDim, TestParamDim);

% Loop through all noisy signals and calculate the deviance between measured parameter value and known true value for both absolute and percentage error
for s = 1:NoisySigDim
    % Calculate the deviance between the signal parameter and the known true value
    AbsParameterDevianceArray(:,s,1) = squeeze(BiexpFitArray(:,s,1))./(squeeze(BiexpFitArray(:,s,1))+squeeze(BiexpFitArray(:,s,2))) - ...
                                    squeeze(BootstrapNFTestParameterArray(:,1))./(squeeze(BootstrapNFTestParameterArray(:,1))+squeeze(BootstrapNFTestParameterArray(:,2))); % SF1 abs error
    AbsParameterDevianceArray(:,s,2) = squeeze(BiexpFitArray(:,s,3)) - squeeze(BootstrapNFTestParameterArray(:,3)); % ADC1 abs error
    AbsParameterDevianceArray(:,s,3) = squeeze(BiexpFitArray(:,s,4)) - squeeze(BootstrapNFTestParameterArray(:,4)); % ADC2 abs error
    PctParameterDevianceArray(:,s,2) = (squeeze(BiexpFitArray(:,s,3)) - squeeze(BootstrapNFTestParameterArray(:,3)))*100./squeeze(BootstrapNFTestParameterArray(:,3)); % ADC1 pct error
    PctParameterDevianceArray(:,s,3) = (squeeze(BiexpFitArray(:,s,4)) - squeeze(BootstrapNFTestParameterArray(:,4)))*100./squeeze(BootstrapNFTestParameterArray(:,4)); % ADC2 pct error
end

% The noise free parameters are already spaced out over the parameter space so no bin sorting is needed
% From the creation file
NumberOfValues = 400;
SF1Values = 0.025:0.05:0.975;
D1D2RatioValues = 2.45:0.9:19.55;
SF1BinDim = length(SF1Values);  D1D2BinDim = length(D1D2RatioValues);  % For a total of 400

NFSF1ValueArray = BootstrapNFTestParameterArray(:,1)./(BootstrapNFTestParameterArray(:,1)+BootstrapNFTestParameterArray(:,2));
NFD1D2RatioArray = BootstrapNFTestParameterArray(:,3)./BootstrapNFTestParameterArray(:,4);

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
% Combine all errors into distribution metrics by parameter
for i = 1:SF1BinDim
    for j = 1:D1D2BinDim
        CurBinIdx = sub2ind([SF1BinDim D1D2BinDim],i,j);
        AbsSF1DevArray = AbsParameterDevianceArray(CurBinIdx,:,1);
        ADC1PctDevArray = PctParameterDevianceArray(CurBinIdx,:,2);
        ADC2PctDevArray = PctParameterDevianceArray(CurBinIdx,:,3);
        SF1ParamEst = SF1ParameterEstimates(CurBinIdx,:);
        ADC1ParamEst = ADC1ParameterEstimates(CurBinIdx,:);
        ADC2ParamEst = ADC2ParameterEstimates(CurBinIdx,:);

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
        
    end
end

PlotFontSize = 18;
ContourFontSize = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SF1, mean absolute error and standard deviation to choose the parameter values below
figure('Position', [50,50,800,800]);
h=pcolor(SF1Values(:), D1D2RatioValues(:), squeeze(BinMeanArray(:,:,1))');set(h,'edgecolor','none');
title('Mean Absolute Error (color) and Std. Dev. (contours) of \itSF_{1}\rm\bf Estimates');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itADC_{1} \rm\bf/ \itADC_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
caxis([-0.3 0.5]);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
colormap(jet); % gives a good contrast for negative numbers
HandleC =  colorbar('location', 'SouthOutside');
% Overlay std dev contours
hold on;
Contours = 0.:0.05:0.4;
% The std dev contours need a bit of averaging for smoothness
flt_size = 3; avg_flt_size = 3; pad_size = flt_size;% window size
AvgFilter = fspecial('average', [avg_flt_size avg_flt_size]);
SF1StdDevPadded = padarray(BinStdDevArray(:,:,1)', [pad_size pad_size], 'replicate'); %zero gradient pad
SF1StdDev = imcrop(filter2(AvgFilter, SF1StdDevPadded), [pad_size + 1, pad_size + 1, SF1BinDim - 1, D1D2BinDim - 1]);
[Cons, HandleSF1Var] = contour(SF1Values(:), D1D2RatioValues(:), SF1StdDev, Contours, 'Linecolor', [.0 .0 .0]);
set(HandleSF1Var, 'ShowText', 'on', 'TextStep', get(HandleSF1Var, 'LevelStep')*200., 'linewidth', 1, 'labelspacing', 300)
HandleText = clabel(Cons,HandleSF1Var);
set(HandleText, 'fontsize', ContourFontSize);   %contour labels font size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select three test cases here to show the bootstrap distributions.  These numbers came from my test set.  Yours will differ.
% Low uncertainty: SF1 = 0.525, D1/D2 = 15.05 (11,15), Index: 291
% Low signal fraction: SF1 = 0.075, D1/D2 = 15.05 (2,15), Index: 282
% Low ratio: SF1 = 0.525, D1/D2 = 2.45  (11,1), Index: 11
% For each noise-free signal show the bootstrap fits for a noisy signal fit close to the true values and one far
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NF Signal Index 291 - Good One
% Noisy Signal 8: 0.5082    0.4859    0.9991    0.0745
% Noisy Signal 13: 0.6690    0.3105    0.6514    0.0355
CurIdx = 291;  NoisySignal1 = 8;  NoisySignal2 = 13;
NoisySignal1Boots = squeeze(BiexpBootFitArray(CurIdx,NoisySignal1,:,:));
NoisySignal2Boots = squeeze(BiexpBootFitArray(CurIdx,NoisySignal2,:,:));
NoisySignal1BiexpFit = squeeze(BiexpFitArray(CurIdx,NoisySignal1,:));
NoisySignal2BiexpFit = squeeze(BiexpFitArray(CurIdx,NoisySignal2,:));
figure('Position', [50,0,1000,800]);
h1 = subplot(4,2,1);histogram(NoisySignal1Boots(:,1),0:0.02:1);title('Noisy Signal 1');ylabel('\itA_{1}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);%xlim([-1 1.5]);
% h1 = subplot(4,2,1);histogram(NoisySignal1Boots(:,1),50);title('Noisy Signal 1');ylabel('A1 Estimates','FontWeight','bold');ylim([0 150]);%xlim([-1 1.5]);
hold on; plot(repmat(0.525,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r'); % Plot true value
hold on; plot(repmat(NoisySignal1BiexpFit(1),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h2 = subplot(4,2,3);histogram(NoisySignal1Boots(:,2),0:0.02:1);ylabel('\itA_{2}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);%xlim([-0.5 2]);
% h2 = subplot(4,2,3);histogram(NoisySignal1Boots(:,2),50);ylabel('A2 Estimates','FontWeight','bold');ylim([0 150]);%xlim([-0.5 2]);
hold on; plot(repmat(0.475,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal1BiexpFit(2),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h3 = subplot(4,2,5);histogram(NoisySignal1Boots(:,3),0:0.08:4);ylabel('\itD_{1}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);%xlim([-2 4]);
% h3 = subplot(4,2,5);histogram(NoisySignal1Boots(:,3),50);ylabel('ADC1 Estimates','FontWeight','bold');ylim([0 150]);%xlim([-2 4]);
hold on; plot(repmat(1,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal1BiexpFit(3),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h4 = subplot(4,2,7);histogram(NoisySignal1Boots(:,4),0:0.004:0.2);ylabel('\itD_{2}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);%xlim([-0.5 0.5]);
% h4 = subplot(4,2,7);histogram(NoisySignal1Boots(:,4),50);ylabel('ADC2 Estimates','FontWeight','bold');ylim([0 150]);%xlim([-0.5 0.5]);
hold on; plot(repmat(0.066,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal1BiexpFit(4),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h5 = subplot(4,2,2);histogram(NoisySignal2Boots(:,1),0:0.02:1);title('Noisy Signal 2');ylim([0 150]);%xlim([-1 1.5]);
% h5 = subplot(4,2,2);histogram(NoisySignal2Boots(:,1),50);title('Noisy Signal 2');ylim([0 150]);%xlim([-1 1.5]);
hold on; plot(repmat(0.525,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal2BiexpFit(1),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h6 = subplot(4,2,4);histogram(NoisySignal2Boots(:,2),0:0.02:1);ylim([0 150]);%xlim([-0.5 2]);
% h6 = subplot(4,2,4);histogram(NoisySignal2Boots(:,2),50);ylim([0 150]);%xlim([-0.5 2]);
hold on; plot(repmat(0.475,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal2BiexpFit(2),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h7 = subplot(4,2,6);histogram(NoisySignal2Boots(:,3),0:0.08:4);ylim([0 300]);%xlim([-2 4]);
% h7 = subplot(4,2,6);histogram(NoisySignal2Boots(:,3),50);ylim([0 300]);%xlim([-2 4]);
hold on; plot(repmat(1,1,2),[0 300], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal2BiexpFit(3),1,2),[0 300], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h8 = subplot(4,2,8);histogram(NoisySignal2Boots(:,4),0:0.004:0.2);ylim([0 150]);%xlim([-0.5 0.5]);
% h8 = subplot(4,2,8);histogram(NoisySignal2Boots(:,4),50);ylim([0 150]);%xlim([-0.5 0.5]);
hold on; plot(repmat(0.066,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal2BiexpFit(4),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NF Signal Index 282 - True Values: A1 = 0.075, A2 = 0.925, D1 = 1, D2 = 0.0664
% Noisy Signal 23: 0.0899    0.9027    1.0071    0.0574
% Noisy Signal 3: 0.5494    0.4398    0.1884    0.0240
squeeze(BiexpFitArray(282,:,:))
CurIdx = 282;  NoisySignal1 = 23;  NoisySignal2 = 3;
NoisySignal1Boots = squeeze(BiexpBootFitArray(CurIdx,NoisySignal1,:,:));
NoisySignal2Boots = squeeze(BiexpBootFitArray(CurIdx,NoisySignal2,:,:));
NoisySignal1BiexpFit = squeeze(BiexpFitArray(CurIdx,NoisySignal1,:));
NoisySignal2BiexpFit = squeeze(BiexpFitArray(CurIdx,NoisySignal2,:));
figure('Position', [50,0,1000,800]);
h1 = subplot(4,2,1);histogram(NoisySignal1Boots(:,1),0:0.02:1);title('Noisy Signal 1');ylabel('\itA_{1}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);%xlim([-1 1.5]);
% h1 = subplot(4,2,1);histogram(NoisySignal1Boots(:,1),50);title('Noisy Signal 1');ylabel('A1 Estimates','FontWeight','bold');ylim([0 150]);%xlim([-1 1.5]);
hold on; plot(repmat(0.075,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r'); % Plot true value
hold on; plot(repmat(NoisySignal1BiexpFit(1),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h2 = subplot(4,2,3);histogram(NoisySignal1Boots(:,2),0:0.02:1);ylabel('\itA_{2}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);%xlim([-0.5 2]);
% h2 = subplot(4,2,3);histogram(NoisySignal1Boots(:,2),50);ylabel('A2 Estimates','FontWeight','bold');ylim([0 150]);%xlim([-0.5 2]);
hold on; plot(repmat(0.925,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal1BiexpFit(2),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h3 = subplot(4,2,5);histogram(NoisySignal1Boots(:,3),0:0.08:4);ylabel('\itD_{1}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);%xlim([-2 4]);
% h3 = subplot(4,2,5);histogram(NoisySignal1Boots(:,3),50);ylabel('ADC1 Estimates','FontWeight','bold');ylim([0 150]);%xlim([-2 4]);
hold on; plot(repmat(1,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal1BiexpFit(3),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h4 = subplot(4,2,7);histogram(NoisySignal1Boots(:,4),0:0.004:0.2);ylabel('\itD_{2}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);%xlim([-0.5 0.5]);
% h4 = subplot(4,2,7);histogram(NoisySignal1Boots(:,4),50);ylabel('ADC2 Estimates','FontWeight','bold');ylim([0 150]);%xlim([-0.5 0.5]);
hold on; plot(repmat(0.066,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal1BiexpFit(4),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h5 = subplot(4,2,2);histogram(NoisySignal2Boots(:,1),0:0.02:1);title('Noisy Signal 2');ylim([0 150]);%xlim([-1 1.5]);
% h5 = subplot(4,2,2);histogram(NoisySignal2Boots(:,1),50);title('Noisy Signal 2');ylim([0 150]);%xlim([-1 1.5]);
hold on; plot(repmat(0.075,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal2BiexpFit(1),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h6 = subplot(4,2,4);histogram(NoisySignal2Boots(:,2),0:0.02:1);ylim([0 150]);%xlim([-0.5 2]);
% h6 = subplot(4,2,4);histogram(NoisySignal2Boots(:,2),50);ylim([0 150]);%xlim([-0.5 2]);
hold on; plot(repmat(0.925,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal2BiexpFit(2),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h7 = subplot(4,2,6);histogram(NoisySignal2Boots(:,3),0:0.08:4);ylim([0 300]);%xlim([-2 4]);
% h7 = subplot(4,2,6);histogram(NoisySignal2Boots(:,3),50);ylim([0 300]);%xlim([-2 4]);
hold on; plot(repmat(1,1,2),[0 300], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal2BiexpFit(3),1,2),[0 300], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h8 = subplot(4,2,8);histogram(NoisySignal2Boots(:,4),0:0.004:0.2);ylim([0 150]);%xlim([-0.5 0.5]);
% h8 = subplot(4,2,8);histogram(NoisySignal2Boots(:,4),50);ylim([0 150]);%xlim([-0.5 0.5]);
hold on; plot(repmat(0.066,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal2BiexpFit(4),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value

% Calculate the 95% CI's
maxBootCI1 = zeros(1,4);
minBootCI1 = zeros(1,4);
for p = 1:4
    curBootCI = prctile(squeeze(NoisySignal1Boots(:,p)),[2.5 97.5]);
    maxBootCI1(p) = curBootCI(2);
    minBootCI1(p) = curBootCI(1);
end

maxBootCI2 = zeros(1,4);
minBootCI2 = zeros(1,4);
for p = 1:4
    curBootCI = prctile(squeeze(NoisySignal2Boots(:,p)),[2.5 97.5]);
    maxBootCI2(p) = curBootCI(2);
    minBootCI2(p) = curBootCI(1);
end

% Overlay the bootstrap CI - MIN
hold on; plot(h1,repmat(minBootCI1(1),1,2),[0 150], 'linewidth',1, 'linestyle','--','Color','b');
hold on; plot(h2,repmat(minBootCI1(2),1,2),[0 150], 'linewidth',1, 'linestyle','--','Color','b');
hold on; plot(h3,repmat(minBootCI1(3),1,2),[0 150], 'linewidth',1, 'linestyle','--','Color','b');
hold on; plot(h4,repmat(minBootCI1(4),1,2),[0 150], 'linewidth',1, 'linestyle','--','Color','b');
hold on; plot(h5,repmat(minBootCI2(1),1,2),[0 150], 'linewidth',1, 'linestyle','--','Color','b');
hold on; plot(h6,repmat(minBootCI2(2),1,2),[0 150], 'linewidth',1, 'linestyle','--','Color','b');
hold on; plot(h7,repmat(minBootCI2(3),1,2),[0 300], 'linewidth',1, 'linestyle','--','Color','b');
hold on; plot(h8,repmat(minBootCI2(4),1,2),[0 150], 'linewidth',1, 'linestyle','--','Color','b');
% And MAX
hold on; plot(h1,repmat(maxBootCI1(1),1,2),[0 150], 'linewidth',1, 'linestyle','--','Color','b');
hold on; plot(h2,repmat(maxBootCI1(2),1,2),[0 150], 'linewidth',1, 'linestyle','--','Color','b');
hold on; plot(h3,repmat(maxBootCI1(3),1,2),[0 150], 'linewidth',1, 'linestyle','--','Color','b');
hold on; plot(h4,repmat(maxBootCI1(4),1,2),[0 150], 'linewidth',1, 'linestyle','--','Color','b');
hold on; plot(h5,repmat(maxBootCI2(1),1,2),[0 150], 'linewidth',1, 'linestyle','--','Color','b');
hold on; plot(h6,repmat(maxBootCI2(2),1,2),[0 150], 'linewidth',1, 'linestyle','--','Color','b');
hold on; plot(h7,repmat(maxBootCI2(3),1,2),[0 300], 'linewidth',1, 'linestyle','--','Color','b');
hold on; plot(h8,repmat(maxBootCI2(4),1,2),[0 150], 'linewidth',1, 'linestyle','--','Color','b');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NF Signal Index 11 - True Values: A1 = 0.525, A2 = 0.475, D1 = 1, D2 = 0.4082
% Noisy Signal 11:  0.5066    0.5092    1.0955    0.3650  % fits are close
% Noisy Signal 15:  0.9842    0.0202    0.6797    0.0000  % fits are NOT close
CurIdx = 11;  NoisySignal1 = 11;  NoisySignal2 = 15;
squeeze(BiexpFitArray(CurIdx,:,:))
NoisySignal1Boots = squeeze(BiexpBootFitArray(CurIdx,NoisySignal1,:,:));
NoisySignal2Boots = squeeze(BiexpBootFitArray(CurIdx,NoisySignal2,:,:));
NoisySignal1BiexpFit = squeeze(BiexpFitArray(CurIdx,NoisySignal1,:));
NoisySignal2BiexpFit = squeeze(BiexpFitArray(CurIdx,NoisySignal2,:));
figure('Position', [0,0,1000,1000]);
h1 = subplot(4,2,1);histogram(NoisySignal1Boots(:,1),0:0.02:1);title('Noisy Signal 1');ylabel('\itA_{1}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);xlim([0 1]);%xlim([-1 1.5]);
% h1 = subplot(4,2,1);histogram(NoisySignal1Boots(:,1),0:0.03:1.5);title('Noisy Signal 1');ylabel('\itA_{1}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);xlim([-1 2]);% CIs
% h1 = subplot(4,2,1);histogram(NoisySignal1Boots(:,1),50);title('Noisy Signal 1');ylabel('A1 Estimates','FontWeight','bold');ylim([0 150]);xlim([-1 1.5]);
hold on; plot(repmat(0.525,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r'); % Plot true value
hold on; plot(repmat(NoisySignal1BiexpFit(1),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h2 = subplot(4,2,3);histogram(NoisySignal1Boots(:,2),0:0.02:1);ylabel('\itA_{2}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);xlim([0 1]);%xlim([-0.5 2]);
% h2 = subplot(4,2,3);histogram(NoisySignal1Boots(:,2),0:0.03:1.5);ylabel('\itA_{2}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);xlim([-1 2]);% CIs
% h2 = subplot(4,2,3);histogram(NoisySignal1Boots(:,2),50);ylabel('A2 Estimates','FontWeight','bold');ylim([0 150]);xlim([-0.5 2]);
hold on; plot(repmat(0.475,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal1BiexpFit(2),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h3 = subplot(4,2,5);histogram(NoisySignal1Boots(:,3),0:0.08:4);ylabel('\itD_{1}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);xlim([0 4]);%xlim([-2 4]);
% h3 = subplot(4,2,5);histogram(NoisySignal1Boots(:,3),0:0.16:8);ylabel('\itD_{1}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);xlim([-2 8]);
% h3 = subplot(4,2,5);histogram(NoisySignal1Boots(:,3),50);ylabel('ADC1 Estimates','FontWeight','bold');ylim([0 150]);xlim([-2 8]);
hold on; plot(repmat(1,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal1BiexpFit(3),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h4 = subplot(4,2,7);histogram(NoisySignal1Boots(:,4),0:0.016:0.8);ylabel('\itD_{2}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);xlim([0 0.8]);%xlim([-0.5 0.5]);
% h4 = subplot(4,2,7);histogram(NoisySignal1Boots(:,4),0:0.016:0.8);ylabel('\itD_{2}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);xlim([-0.4 1]);%CIs
% h4 = subplot(4,2,7);histogram(NoisySignal1Boots(:,4),50);ylabel('ADC2 Estimates','FontWeight','bold');ylim([0 150]);xlim([-0.5 0.5]);
hold on; plot(repmat(0.4082,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal1BiexpFit(4),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h5 = subplot(4,2,2);histogram(NoisySignal2Boots(:,1),0:0.02:1);title('Noisy Signal 2');ylim([0 150]);xlim([0 1]);%xlim([-1 1.5]);
% h5 = subplot(4,2,2);histogram(NoisySignal2Boots(:,1),0:0.03:1.5);title('Noisy Signal 2');ylim([0 150]);xlim([0 1.5]);%CIs
% h5 = subplot(4,2,2);histogram(NoisySignal2Boots(:,1),50);title('Noisy Signal 2');ylim([0 150]);xlim([-1 1.5]);
hold on; plot(repmat(0.525,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal2BiexpFit(1),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h6 = subplot(4,2,4);histogram(NoisySignal2Boots(:,2),0:0.02:1);ylim([0 150]);xlim([0 1]);%xlim([-0.5 2]);
% h6 = subplot(4,2,4);histogram(NoisySignal2Boots(:,2),0:0.03:1.5);ylim([0 150]);xlim([-0.5 1.5]);%CIs
% h6 = subplot(4,2,4);histogram(NoisySignal2Boots(:,2),50);ylim([0 150]);xlim([-0.5 2]);
hold on; plot(repmat(0.475,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal2BiexpFit(2),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h7 = subplot(4,2,6);histogram(NoisySignal2Boots(:,3),0:0.08:4);ylim([0 300]);xlim([0 4]);
% h7 = subplot(4,2,6);histogram(NoisySignal2Boots(:,3),0:0.16:8);ylim([0 300]);xlim([0 8]);%CIs
% h7 = subplot(4,2,6);histogram(NoisySignal2Boots(:,3),50);ylim([0 300]);xlim([-2 8]);
hold on; plot(repmat(1,1,2),[0 300], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal2BiexpFit(3),1,2),[0 300], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value
h8 = subplot(4,2,8);histogram(NoisySignal2Boots(:,4),0:0.016:0.8);ylim([0 150]);xlim([0 0.8]);
% h8 = subplot(4,2,8);histogram(NoisySignal2Boots(:,4),0:0.016:0.8);ylim([0 150]);xlim([-0.2 0.8]);%CIs
% h8 = subplot(4,2,8);histogram(NoisySignal2Boots(:,4),50);ylim([0 150]);xlim([-0.5 0.5]);
hold on; plot(repmat(0.4082,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','r');
hold on; plot(repmat(NoisySignal2BiexpFit(4),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g'); % Plot original fit value

% Calculate the 95% CI's
minBootCI1 = zeros(1,4);
maxBootCI1 = zeros(1,4);
for p = 1:4
    curBootCI = prctile(squeeze(NoisySignal1Boots(:,p)),[2.5 97.5]);
    maxBootCI1(p) = curBootCI(2);
    minBootCI1(p) = curBootCI(1);
end

maxBootCI2 = zeros(1,4);
minBootCI2 = zeros(1,4);
for p = 1:4
    curBootCI = prctile(squeeze(NoisySignal2Boots(:,p)),[2.5 97.5]);
    maxBootCI2(p) = curBootCI(2);
    minBootCI2(p) = curBootCI(1);
end

length(find(NoisySignal1Boots(:,3) > 4))

% Overlay the bootstrap CI - MIN
hold on; plot(h1,repmat(minBootCI1(1),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','b');
hold on; plot(h2,repmat(minBootCI1(2),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','b');
hold on; plot(h3,repmat(minBootCI1(3),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','b');
hold on; plot(h4,repmat(minBootCI1(4),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','b');
hold on; plot(h5,repmat(minBootCI2(1),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','b');
hold on; plot(h6,repmat(minBootCI2(2),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','b');
hold on; plot(h7,repmat(minBootCI2(3),1,2),[0 300], 'linewidth',2, 'linestyle','--','Color','b');
hold on; plot(h8,repmat(minBootCI2(4),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','b');
% And MAX
hold on; plot(h1,repmat(maxBootCI1(1),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','b');
hold on; plot(h2,repmat(maxBootCI1(2),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','b');
hold on; plot(h3,repmat(maxBootCI1(3),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','b');
hold on; plot(h4,repmat(maxBootCI1(4),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','b');
hold on; plot(h5,repmat(maxBootCI2(1),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','b');
hold on; plot(h6,repmat(maxBootCI2(2),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','b');
hold on; plot(h7,repmat(maxBootCI2(3),1,2),[0 300], 'linewidth',2, 'linestyle','--','Color','b');
hold on; plot(h8,repmat(maxBootCI2(4),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','b');

% What are the estimate confidence intervals from the regression fit? Compare to t-dist estimates
Jacobian1 = squeeze(BiexpFitJacobianArray(11,11,:,:));
Jacobian2 = squeeze(BiexpFitJacobianArray(11,15,:,:));
Residuals1 = squeeze(BiexpResidualArray(11,11,:));
Residuals2 = squeeze(BiexpResidualArray(11,15,:));
BDim = 11;
DOFModel = BDim-ParamDim;
alpha = 0.05; % 95% CI
TFactor = tinv(1-alpha/2,DOFModel);

% Signal 1 - convert to parameter standard errors
[Q1,R1] = qr(Jacobian1,0);
Rinv1 = R1 \ eye(size(R1));
JTJInv1 = Rinv1*Rinv1';
% Get standard error value
RSS1 = sum(Residuals1.^2);
StdVariance1 = RSS1/DOFModel;
CovMx1 = StdVariance1*JTJInv1;
ParamSE1 = sqrt(diag(CovMx1));
MinCI1 = NoisySignal1BiexpFit - ParamSE1*TFactor;
MaxCI1 = NoisySignal1BiexpFit + ParamSE1*TFactor;

% Signal 2 - convert to parameter standard errors
[Q2,R2] = qr(Jacobian2,0);
Rinv2 = R2 \ eye(size(R2));
JTJInv2 = Rinv2*Rinv2';
% Get standard error value
RSS2 = sum(Residuals2.^2);
StdVariance2 = RSS2/DOFModel;
CovMx2 = StdVariance2*JTJInv2;
ParamSE2 = sqrt(diag(CovMx2));
MinCI2 = NoisySignal2BiexpFit - ParamSE2*TFactor;
MaxCI2 = NoisySignal2BiexpFit + ParamSE2*TFactor;

% Overlay the t-distributions - MIN
hold on; plot(h1,repmat(MinCI1(1),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','c');
hold on; plot(h2,repmat(MinCI1(2),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','c');
hold on; plot(h3,repmat(MinCI1(3),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','c');
hold on; plot(h4,repmat(MinCI1(4),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','c');
hold on; plot(h5,repmat(MinCI2(1),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','c');
hold on; plot(h6,repmat(MinCI2(2),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','c');
hold on; plot(h7,repmat(MinCI2(3),1,2),[0 300], 'linewidth',2, 'linestyle','--','Color','c');
hold on; plot(h8,repmat(MinCI2(4),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','c');
% And MAX
hold on; plot(h1,repmat(MaxCI1(1),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','c');
hold on; plot(h2,repmat(MaxCI1(2),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','c');
hold on; plot(h3,repmat(MaxCI1(3),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','c');
hold on; plot(h4,repmat(MaxCI1(4),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','c');
hold on; plot(h5,repmat(MaxCI2(1),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','c');
hold on; plot(h6,repmat(MaxCI2(2),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','c');
hold on; plot(h7,repmat(MaxCI2(3),1,2),[0 300], 'linewidth',2, 'linestyle','--','Color','c');
hold on; plot(h8,repmat(MaxCI2(4),1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','c');
