% Chap3NLLSKurtosisFitToMonoexpData.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap3NLLSKurtosisFitToMonoexpData()

load('YourPath\MonoexpNoisySignals_SNR25_11BValues.mat'); % HPC Scratch
[NFSigDim NoisySigDim BDim] = size(NoisySignalArray);

% Add BValueArray back to itself, as cluster parfor can't see it - MATLAB strangeness
TmpBValues = BValueArray;
BValueArray = TmpBValues;

% Calculate the noise std. dev.
TrueSTDNoise = 0.04 % For SNR of 25 = 0.04

% Create the array to store the fitted parameters based on the three spatial dimensions plus the axes and parameters.
% Make these all 1D for the parfor array and then reshape after.
KurtFitArray = zeros(NFSigDim, NoisySigDim, 3); % Amp, ADC, K
KurtResidualArray = zeros(NFSigDim, NoisySigDim, BDim);
KurtFitJacobianArray = zeros(NFSigDim, NoisySigDim,BDim, 3);
KurtFitOutputFlagArray = zeros(NFSigDim, NoisySigDim, 1); % Saves the algorithm return code

% Values for multiple start points
MinAmpValue = 0.8;  MaxAmpValue = 1.2;
MinADCValue = 0.5;   MaxADCValue = 1.5;
MinKValue = -0.2;   MaxKValue = 0.2;
ParamDim = 3;
NumberRandomStarts = 5;

FitOpt = optimset('Display', 'off', 'Algorithm', 'trust-region-reflective', 'Jacobian', 'on');

% Now loop through and calculate all values for each voxel.  Flatten the 2D signal array to 1D.
for i = 1:NFSigDim
    parfor j = 1:NoisySigDim
        % FITTING
        VoxelSignal = squeeze(NoisySignalArray(i,j,:))';
        % Remove bias using Gudbjarsson equation
        % For each signal add a random amount of small error to take into account the fact that the std of the noise would not be known
        STDNoise = TrueSTDNoise + randn(1)*TrueSTDNoise*0.1;
        CorrVoxelSignal = sqrt(abs(VoxelSignal.^2 - STDNoise.^2));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % KURTOSIS FIT
        % Setup multi-start points
        FitJac = zeros(NumberRandomStarts,BDim,ParamDim); FitResiduals = zeros(NumberRandomStarts,BDim); FitFlag = zeros(NumberRandomStarts,1);
        FitParams = zeros(NumberRandomStarts,ParamDim); FitRSS = zeros(NumberRandomStarts,1);
        FitRSS(:) = 1e30; % Keep this high, since we are looking for a minimum;
        LB = [0, 0, -1.]; %[0., 0., 0., 0.]; % AMP, ADC, K
        UB = [2*max(CorrVoxelSignal), 4., 2.]; % AMP, ADC, K

        for k = 1:NumberRandomStarts
            StartPoints = zeros(1,3);
            AmpValue = (MaxAmpValue-MinAmpValue)*CorrVoxelSignal(1)*rand(1, 1) + MinAmpValue;
            ADCValue = (MaxADCValue-MinADCValue)*rand(1, 1) + MinADCValue;
            KValue = (MaxKValue-MinKValue)*rand(1, 1) + MinKValue;
            StartPoints(:,1) = AmpValue; StartPoints(:,2) = ADCValue; StartPoints(:,3) = KValue; 

            [RetParams, RSS, Residuals,XFlag,Output,RetLamb,RetJac]  =  lsqcurvefit(@KurtDecayWithJac, StartPoints, BValueArray, CorrVoxelSignal, LB, UB, FitOpt);

            FitParams(k,:) = [RetParams(1), RetParams(2), RetParams(3)];
            FitResiduals(k,:) = Residuals;  FitRSS(k) = RSS;  FitFlag(k) = XFlag;
            FitJac(k,:,:) = full(RetJac);
        end

        % Go through and find the minimum RSS value and set those values in the array
        MinRSSIndex = find(FitRSS == min(FitRSS), 1 ); % wrap with min in case the find returns multiple values
        KurtFitArray(i,j,:) = FitParams(MinRSSIndex,:); % Amp, ADC, K
        KurtResidualArray(i,j,:) = FitResiduals(MinRSSIndex,:); % Based on BDim
        KurtFitJacobianArray(i,j,:,:) = FitJac(MinRSSIndex,:,:); % Must be 3x3 to accommodate parameter Jacobians
        KurtFitOutputFlagArray(i,j) = FitFlag(MinRSSIndex); % Saves the algorithm return code
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    i
    datestr(now)
end

% Save all data in separate values
% Large files must be saved with the '-v7.3' flag
save('YourPath\KurtFitArray11BValuesKurtChap.mat','KurtFitArray','-v7.3');
save('YourPath\KurtResidualArray11BValuesKurtChap.mat','KurtResidualArray','-v7.3');
save('YourPath\KurtFitJacobianArray11BValuesKurtChap.mat','KurtFitJacobianArray','-v7.3');
save('YourPath\KurtFitOutputFlagArray11BValuesKurtChap.mat','KurtFitOutputFlagArray','-v7.3');

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