% Chap3DisplayBootstrapHistograms.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap3DisplayBootstrapHistograms()

% This function for the paper goes through and adds noise for the noisy signals to be used in the Monte Carlo simulation.
% It also creates the noise free signals based on the particular range of b values for the noise addition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the appropriate signal
load('YourPath\BootstrapNoiseFreeParameters.mat');

% Test signals SF1 = 0.025, D1/D2 ratio = 15.05
SF1Test = 1; D1D2RatioTest = 15;
TestSignalIdx = sub2ind([20 20],SF1Test,D1D2RatioTest); % Should be 281

load('YourPath\BootstrapNoisySignals_SNR25.mat'); % For local testing
% Expecting this to have BValueArray and BootstrapNoisySignalArray
[NFSigDim NoisySigDim BDim] = size(BootstrapNoisySignalArray);
TestBootstrapSignalArray = squeeze(BootstrapNoisySignalArray(TestSignalIdx,:,:));

% From this set, a noisy fit with little kurtosis is selected (since the true value is close to monoexponential)
SelectedNoisyFitIdx = 20; % Amp: 0.9672  Dapp: 0.0712  Kapp: 0.0681

% Calculate the noise std. dev.
TrueSTDNoise = 0.04 % For SNR of 25 = 0.04

% Number of bootstrap samples per measurement
BootstrapSamplesPerMeasurement = 1000;  BootDim = BootstrapSamplesPerMeasurement;

% Declare fit result arrays
KurtFitArray = zeros(1, 3); % Amp, ADC, K
KurtResidualArray = zeros(1, BDim);
KurtFitJacobianArray = zeros(BDim, 3);
KurtFitOutputFlagArray = zeros(1, 1); % Saves the algorithm return code

% Extra bootstrap info
KurtBootFitArray = zeros(BootDim, 3); % Amp, ADC, K
KurtBootResidualArray = zeros(BootDim, BDim);
KurtBootFitJacobianArray = zeros(BootDim, BDim, 3);
KurtBootFitOutputFlagArray = zeros(BootDim, 1);

FitOpt = optimset('Display', 'off', 'Algorithm', 'trust-region-reflective', 'Jacobian', 'on');

% Values for multiple start points
MinAmpValue = 0.8;  MaxAmpValue = 1.2;
MinADCValue = 0.5;   MaxADCValue = 1.5;
MinKValue = -0.2;   MaxKValue = 0.2;
ParamDim = 3;
NumberRandomStarts = 5;

% FITTING with 2SD
VoxelSignal = squeeze(TestBootstrapSignalArray(SelectedNoisyFitIdx,:));

% Remove bias using Gudbjarsson equation
% For each signal add a random amount of small error to take into account the fact that the std of the noise would not be known
STDNoise = TrueSTDNoise + randn(1)*TrueSTDNoise*0.1;
CorrVoxelSignal = sqrt(abs(VoxelSignal.^2 - STDNoise.^2));
MaxSignalIndex = BDim;
% Use this max signal index to limit the b values and signal
CurCorrVoxelSignal = CorrVoxelSignal(1:MaxSignalIndex);
CurBValueArray = BValueArray(1:MaxSignalIndex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KURTOSIS FIT
% Setup multi-start points
FitJac = zeros(NumberRandomStarts,MaxSignalIndex,ParamDim); FitResiduals = zeros(NumberRandomStarts,MaxSignalIndex); FitFlag = zeros(NumberRandomStarts,1);
FitParams = zeros(NumberRandomStarts,ParamDim); FitRSS = zeros(NumberRandomStarts,1);
FitRSS(:) = 1e30; % Keep this high, since we are looking for a minimum;
LB = [0, 0, -1.]; %[0., 0., 0., 0.]; % AMP, ADC, K
UB = [2*max(CurCorrVoxelSignal), 4., 2.]; % AMP, ADC, K

for k = 1:NumberRandomStarts
    StartPoints = zeros(1,3);
    AmpValue = (MaxAmpValue-MinAmpValue)*CurCorrVoxelSignal(1)*rand(1, 1) + MinAmpValue;
    ADCValue = (MaxADCValue-MinADCValue)*rand(1, 1) + MinADCValue;
    KValue = (MaxKValue-MinKValue)*rand(1, 1) + MinKValue;
    StartPoints(:,1) = AmpValue; StartPoints(:,2) = ADCValue; StartPoints(:,3) = KValue; 

    [RetParams, RSS, Residuals,XFlag,Output,RetLamb,RetJac]  =  lsqcurvefit(@KurtDecayWithJac, StartPoints, CurBValueArray, CurCorrVoxelSignal, LB, UB, FitOpt);

    FitParams(k,:) = [RetParams(1), RetParams(2), RetParams(3)];
    FitResiduals(k,:) = Residuals;  FitRSS(k) = RSS;  FitFlag(k) = XFlag;
    FitJac(k,:,:) = full(RetJac);
end

% Go through and find the minimum RSS value and set those values in the array
MinRSSIndex = find(FitRSS == min(FitRSS), 1 ); % wrap with min in case the find returns multiple values

KurtFitArray(:) = FitParams(MinRSSIndex,:); % Amp, ADC, K
KurtResidualArray(:) = FitResiduals(MinRSSIndex,:); % Based on BDim
KurtFitJacobianArray(:,:) = FitJac(MinRSSIndex,:,:); % Must be 3x3 to accommodate parameter Jacobians
KurtFitOutputFlagArray(:) = FitFlag(MinRSSIndex); % Saves the algorithm return code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup the bootstrap samples for this biexponential fit.  Resample the residuals for parametric bootstrap
% Help from here: http://www.math.ubc.ca/~keshet/MCB2012/SlidesDodo/DataFitLect3.pdf
CurResiduals = FitResiduals(MinRSSIndex,:);
% Find the y values that correspond to this fit
amp = FitParams(MinRSSIndex,1); adc = FitParams(MinRSSIndex,2); k = FitParams(MinRSSIndex,3);
InlineBiexp = @(xdata) amp.*exp(-adc.*xdata + (k.*adc.^2.*xdata.^2)/6);
yFitValues = InlineBiexp(BValueArray);
% Setup the random residual samples with replacement
% Need twice as many since two measurements will be needed for the Rician bias
[~, bootIndices] = bootstrp(2.*BootstrapSamplesPerMeasurement, [], CurResiduals);
bootResiduals = CurResiduals(bootIndices)'; % reshape for addition
% Generate data sets based on the sampled residuals
% Break out the residuals into two groups for Rician noise
FirstBootResiduals = bootResiduals(1:BootstrapSamplesPerMeasurement,:); 
SecondBootResiduals = bootResiduals(BootstrapSamplesPerMeasurement+1:2*BootstrapSamplesPerMeasurement,:);
ySignals = repmat(yFitValues, BootstrapSamplesPerMeasurement, 1);
% Now use Rician calculations to make Gaussian residual perturbations Rician
yBoot = sqrt((ySignals + FirstBootResiduals).^2 + (SecondBootResiduals).^2);
% With the bootstrap residuals randomly selected, loop through and fit each sample
for k = 1:BootDim
    % Retrive the bootstrap sample data
    CurData = yBoot(k,:);
    StartPoints = [1., 1, 0.];
    LB = [0., 0., -1.]; UB = [2*max(CurCorrVoxelSignal), 4, 2];
    [RetParams, ~, Residuals, FitFlag, ~, ~, RetJac]  =  lsqcurvefit(@KurtDecayWithJac, StartPoints, BValueArray, CurData, LB, UB, FitOpt);
    KurtBootFitArray(k,:) = [RetParams(1), RetParams(2), RetParams(3)];
    KurtBootResidualArray(k,:) = Residuals;
    KurtBootFitJacobianArray(k,:,:) = full(RetJac);
    KurtBootFitOutputFlagArray(k) = FitFlag;
    k
end

S0 = 0.97, Dapp = 0.071, Kapp = 0.068

figure('Position', [0,0,1000,800]);
h1 = subplot(2,2,1);histogram(KurtBootFitArray(:,1),50);ylabel('\itS_{0}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);%xlim([-1 1.5]);
hold on; plot(repmat(0.97,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g');
h2 = subplot(2,2,2);histogram(KurtBootFitArray(:,2),50);ylabel('\itD_{app}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);%xlim([-1 1.5]);
hold on; plot(repmat(0.071,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g');
h3 = subplot(2,2,3);histogram(KurtBootFitArray(:,3),50);ylabel('\itK_{app}\rm\bf Estimates','FontWeight','bold');ylim([0 150]);%xlim([-1 1.5]);
hold on; plot(repmat(00.068,1,2),[0 150], 'linewidth',2, 'linestyle','--','Color','g');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signal,jacobian] = KurtDecayWithJac(x, xdata)
amp = x(1); adc = x(2); k = x(3);
signal = amp.*exp(-adc.*xdata + (k.*adc.^2.*xdata.^2)/6);
% jacobian = (1/3.*amp.*adc).*(adc*k.*xdata-3).*exp(1/6.*adc.*xdata.*(adc*k.*xdata-6)); % This is d/dx
JacobPt1 = exp(1/6.*adc.^2*k.*xdata.^2-adc.*xdata); % This is d/d(amp)
JacobPt2 = (1/3.*amp.*xdata).*(adc*k.*xdata-3).*exp(1/6.*adc.*xdata.*(adc*k.*xdata-6)); % This is d/d(adc)
JacobPt3 = (1/6.*amp.*adc.^2.*xdata.^2).*exp(1/6.*adc.*xdata.*(adc*k.*xdata-6)); % This is d/k
jacobian = [JacobPt1' JacobPt2' JacobPt3']; % Complete Jacobian       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%