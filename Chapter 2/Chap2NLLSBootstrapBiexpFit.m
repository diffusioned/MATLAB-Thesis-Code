% Chap2NLLSBootstrapBiexpFit.m
% MATLAB file for fitting simulated bootstrap data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap2NLLSBootstrapBiexpFit()

load('YourPath\BootstrapNoisySignals_SNR25.mat'); % For local testing
% Expecting this to have BValueArray and BootstrapNoisySignalArray
[NFSigDim NoisySigDim BDim] = size(BootstrapNoisySignalArray);

% Calculate the noise std. dev.
TrueSTDNoise = 0.04 % For SNR of 25 = 0.04

% Number of bootstrap samples per measurement
BootstrapSamplesPerMeasurement = 1000;  BootDim = BootstrapSamplesPerMeasurement;

% Declare fit result arrays
BiexpFitArray = zeros(NFSigDim, NoisySigDim, 4); % A1,A2,D1,D2
BiexpResidualArray = zeros(NFSigDim, NoisySigDim, BDim);
BiexpFitJacobianArray = zeros(NFSigDim, NoisySigDim, BDim,4);
BiexpFitOutputFlagArray = zeros(NFSigDim, NoisySigDim, 1); % Saves the algorithm return code
MonoexpFitArray = zeros(NFSigDim, NoisySigDim, 2); % A0,ADC
MonoexpResidualArray = zeros(NFSigDim, NoisySigDim, BDim);
% Extra bootstrap info
BiexpBootFitArray = zeros(NFSigDim, NoisySigDim, BootDim, 4); % A1,A2,D1,D2
BiexpBootResidualArray = zeros(NFSigDim, NoisySigDim, BootDim, BDim);
BiexpBootFitJacobianArray = zeros(NFSigDim, NoisySigDim, BootDim, BDim,4);
BiexpBootFitOutputFlagArray = zeros(NFSigDim, NoisySigDim, BootDim, 1);

% Values for multiple start points
MinSF1Value = 0.0;  MaxSF1Value = 1;
MinADCValue = 0.001;   MaxADCValue = 4;
ParamDim = 4;
NumberRandomStarts = 5;

FitOpt = optimset('Display', 'off', 'Algorithm', 'trust-region-reflective', 'Jacobian', 'on');

% To use the parfor array, the BValueArray has to be made local.  Unknown file access bug?  MATLAB magic?
TmpBValueArray = BValueArray;
BValueArray = TmpBValueArray;

% Fit the signal to the original data value
for i = 1:NFSigDim
    for j = 1: NoisySigDim
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        VoxelSignal = squeeze(BootstrapNoisySignalArray(i,j,:))';
        % Remove bias using Gudbjarsson equation
        % For each signal add a random amount of small error to take into account the fact that the std of the noise would not be known
        STDNoise = TrueSTDNoise + randn(1)*TrueSTDNoise*0.1;
        CorrVoxelSignal = sqrt(abs(VoxelSignal.^2 - STDNoise.^2));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MONOEXPONENTIAL FIT
        % Construct the starting guess and bounding vectors
        StartGuess = [max(CorrVoxelSignal), 1.5];% 2 fits
        LB = [0., 0.];    UB = [2*max(CorrVoxelSignal), 4.];
        [FitParams, ~, Residuals]  =  lsqcurvefit(@MonoexpDecayWithJac, StartGuess, BValueArray, CorrVoxelSignal, LB, UB, FitOpt);
        MonoexpFitArray(i,j,:) = [FitParams(1), FitParams(2)];
        MonoexpResidualArray(i,j,:) = Residuals;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BIEXPONENTIAL FIT
        % Only bootstrapping and multi-starting for the biexponential, since it deviates more
        FitJac = zeros(NumberRandomStarts,BDim,ParamDim); FitResiduals = zeros(NumberRandomStarts,BDim); FitFlag = zeros(NumberRandomStarts,1);
        FitParams = zeros(NumberRandomStarts,ParamDim); FitRSS = zeros(NumberRandomStarts,1);
        FitRSS(:) = 1e30; % Keep this high, since we are looking for a minimum;
        LB = [0, 0, 0., 0.]; %[0., 0., 0., 0.]; % with constraints on linear
        UB = [max(CorrVoxelSignal), max(CorrVoxelSignal), 4., 4.]; % [GuessA0, 0.004, GuessA0, 0.004]; % with constraints on linear

        for k = 1:NumberRandomStarts
            StartPoints = zeros(1,4);
            SignalFraction1Value = (MaxSF1Value-MinSF1Value).*rand(1, 1) + MinSF1Value; % SF1 = A1 here because TrueD1Value = 1
            Amp1Value = SignalFraction1Value*CorrVoxelSignal(1);
            Amp2Value = (1 - SignalFraction1Value) * CorrVoxelSignal(1);
            D1Value = (MaxADCValue-MinADCValue).*rand(1, 1) + MinADCValue;
            D2Value = (MaxADCValue-MinADCValue).*rand(1, 1) + MinADCValue;
            StartPoints(:,1) = Amp1Value; StartPoints(:,2) = Amp2Value;
            StartPoints(:,3) = D1Value; StartPoints(:,4) = D2Value;
            [RetParams, RSS, Residuals, XFlag, ~, ~, RetJac]  =  lsqcurvefit(@BiexpDecayIndAmpsWithJac, StartPoints, BValueArray, CorrVoxelSignal, LB, UB, FitOpt);
            if RetParams(3) < RetParams(4)
                % D2 is actually the higher diffusion value
                FitParams(k,:) = [RetParams(2), RetParams(1), RetParams(4), RetParams(3)];
            else
                FitParams(k,:) = [RetParams(1), RetParams(2), RetParams(3), RetParams(4)];
            end
            FitResiduals(k,:) = Residuals;  FitRSS(k) = RSS;  FitFlag(k) = XFlag;
            FitJac(k,:,:) = full(RetJac);
        end
        % Go through and find the minimum RSS value and set those values in the array
        MinRSSIndex = find(FitRSS == min(FitRSS), 1 ); % wrap with min in case the find returns multiple values
        BiexpFitArray(i,j,:) = FitParams(MinRSSIndex,:); % A1,A2,D1,D2
        BiexpResidualArray(i,j,:) = FitResiduals(MinRSSIndex,:); % Based on BDim
        BiexpFitJacobianArray(i,j,:,:) = FitJac(MinRSSIndex,:,:); % Must be 4x4 to accommodate parameter Jacobians
        BiexpFitOutputFlagArray(i,j) = FitFlag(MinRSSIndex); % Saves the algorithm return code
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup the bootstrap samples for this biexponential fit.  Resample the residuals for parametric bootstrap
        % Help from here: http://www.math.ubc.ca/~keshet/MCB2012/SlidesDodo/DataFitLect3.pdf
        CurResiduals = FitResiduals(MinRSSIndex,:);
        % Find the y values that correspond to this fit
        a1 = FitParams(MinRSSIndex,1); a2 = FitParams(MinRSSIndex,2); adc1 = FitParams(MinRSSIndex,3);  adc2 = FitParams(MinRSSIndex,4);
        InlineBiexp = @(xdata) a1*exp(-adc1*xdata) + a2*exp(-adc2*xdata);
        yFitValues = InlineBiexp(BValueArray);
        % Setup the random residual samples with replacement
        [~, bootIndices] = bootstrp(BootstrapSamplesPerMeasurement, [], CurResiduals);
        bootResiduals = CurResiduals(bootIndices);
        % Generate data sets based on the sampled residuals
%         yBoot = zeros(length(BValueArray), BootstrapSamplesPerMeasurement);
        bootResiduals = bootResiduals'; % reshape for addition
        yBoot = repmat(yFitValues, BootstrapSamplesPerMeasurement, 1) + bootResiduals;
        % With the bootstrap residuals randomly selected, loop through and fit each sample
        parfor k = 1:BootDim
            % Retrive the bootstrap sample data
            CurData = yBoot(k,:);
            StartPoints = [0.5, 0.5, 1, 1/6.];
            LB = [0., 0., 0., 0.]; UB = [Inf, Inf, Inf, Inf];
            [RetParams, ~, Residuals, FitFlag, ~, ~, RetJac]  =  lsqcurvefit(@BiexpDecayIndAmpsWithJac, StartPoints, BValueArray, CurData, LB, UB, FitOpt);
            if RetParams(3) < RetParams(4)
                % D2 is actually the higher diffusion value
                BiexpBootFitArray(i,j,k,:) = [RetParams(2), RetParams(1), RetParams(4), RetParams(3)];
            else
                BiexpBootFitArray(i,j,k,:) = [RetParams(1), RetParams(2), RetParams(3), RetParams(4)];
            end
            BiexpBootResidualArray(i,j,k,:) = Residuals;
            BiexpBootFitJacobianArray(i,j,k,:,:) = full(RetJac);
            BiexpBootFitOutputFlagArray(i,j,k) = FitFlag;
        end
    end
    i
    datestr(now)
end

% Save all data in separate values
% Large files must be saved with the '-v7.3' flag
save('YourPath\BiexpFitArray.mat','BiexpFitArray','-v7.3');
save('YourPath\BiexpResidualArray.mat','BiexpResidualArray','-v7.3');
save('YourPath\BiexpFitJacobianArray.mat','BiexpFitJacobianArray','-v7.3');
save('YourPath\BiexpFitOutputFlagArray.mat','BiexpFitOutputFlagArray','-v7.3');
save('YourPath\MonoexpFitArray.mat','MonoexpFitArray','-v7.3');
save('YourPath\MonoexpResidualArray.mat','MonoexpResidualArray','-v7.3');
save('YourPath\BiexpBootFitArray.mat','BiexpBootFitArray','-v7.3');
save('YourPath\BiexpBootResidualArray.mat','BiexpBootResidualArray','-v7.3');
save('YourPath\BiexpBootFitJacobianArray.mat','BiexpBootFitJacobianArray','-v7.3');
save('YourPath\BiexpBootFitOutputFlagArray.mat','BiexpBootFitOutputFlagArray','-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signal,jacobian] = MonoexpDecayWithJac(x, xdata)
amp  = x(1);adc = x(2);
signal = amp*exp(-adc*xdata);
                    
% Need to create a Jacobian from the current data
JacobPt1 = exp(-adc.*xdata);
JacobPt2 = -amp.*xdata.*exp(-adc.*xdata);
jacobian = [JacobPt1' JacobPt2'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signal,jacobian] = BiexpDecayIndAmpsWithJac(params, xdata)
a1  = params(1); a2 = params(2); adc1 = params(3); adc2 = params(4);
signal = a1*exp(-adc1*xdata) + a2*exp(-adc2*xdata);

% Need to create a Jacobian from the current data
JacobPt1 = exp(-adc1.*xdata);
JacobPt2 = exp(-adc2.*xdata);
JacobPt3 = -a1.*xdata.*exp(-adc1.*xdata);
JacobPt4 = -a2.*xdata.*exp(-adc2.*xdata);
jacobian = [JacobPt1' JacobPt2' JacobPt3' JacobPt4'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%