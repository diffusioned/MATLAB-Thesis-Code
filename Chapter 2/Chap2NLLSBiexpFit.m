% Chap2NLLSBiexpFit.m
% MATLAB file for fitting simulated data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap2NLLSBiexpFit()
load('YourPath\YourNoisySignals_SNR25.mat'); % Your noisy signal data
[NFSigDim NoisySigDim BDim] = size(NoisySignalArray); % Get the dimensions of the noisy data array
% Your noisy file above should contain the NoisySignalArray and the BValueArray.  Both are needed.

TrueSTDNoise = 0.04 % For SNR of 25 = 0.04

% Create the array to store the fitted parameters based on the three spatial dimensions plus the axes and parameters.
% Make these all 1D for the parfor array and then reshape after.
BiexpFitArray = zeros(NFSigDim, NoisySigDim,4); % A1,A2,D1,D2
BiexpResidualArray = zeros(NFSigDim, NoisySigDim, BDim);
BiexpFitJacobianArray = zeros(NFSigDim, NoisySigDim,BDim,4); % Must be 4x4 to accommodate parameter Jacobians
BiexpFitOutputFlagArray = zeros(NFSigDim, NoisySigDim,1); % Saves the algorithm return code
MonoexpFitArray = zeros(NFSigDim, NoisySigDim, 3); % A0,ADC,RSS
MonoexpResidualArray = zeros(NFSigDim, NoisySigDim, BDim);

% Values for multiple start points
MinSF1Value = 0.0;  MaxSF1Value = 1;
MinADCValue = 0.001;   MaxADCValue = 4;
ParamDim = 4;
NumberRandomStarts = 5;

FitOpt = optimset('Display', 'off', 'Algorithm', 'trust-region-reflective', 'Jacobian', 'on');

% Now loop through and calculate all values for each voxel.  Flatten the 2D signal array to 1D.
for i = 1:NFSigDim
    parfor j = 1:NoisySigDim
        % FITTING
        VoxelSignal = squeeze(NoisySignalArray(i,j,:))';

        % Remove bias using Gudbjarsson and Patz equation from "The Rician Distribution of Noisy Data", 1995
        % For each signal add a random amount of small error to take into account the fact that the std of the noise would not be known
        STDNoise = TrueSTDNoise + randn(1)*TrueSTDNoise*0.1;
        CorrVoxelSignal = sqrt(abs(VoxelSignal.^2 - STDNoise.^2));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MONOEXPONENTIAL FIT
        % Construct the starting guess and bounding vectors
        StartGuess = [max(CorrVoxelSignal), 1.5];% 2 fits
        LB = [0., 0.];    UB = [2*max(CorrVoxelSignal), 4.];
        [FitParams, RSS, Residuals]  =  lsqcurvefit(@MonoexpDecayWithJac, StartGuess, BValueArray, CorrVoxelSignal, LB, UB, FitOpt);
        MonoexpFitArray(i,j,:) = [FitParams(1), FitParams(2), RSS];
        MonoexpResidualArray(i,j,:) = Residuals;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BIEXPONENTIAL FIT
        % Setup multi-start points
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

            [RetParams, RSS, Residuals,XFlag,Output,RetLamb,RetJac]  =  lsqcurvefit(@BiexpDecayIndAmpsWithJac, StartPoints, BValueArray, CorrVoxelSignal, LB, UB, FitOpt);

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
    end
    i
    datestr(now)
end

% Save all data in separate .mat files
% Large files must be saved with the '-v7.3' flag
save('YourPath\BiexpFitArray.mat','BiexpFitArray','-v7.3');
save('YourPath\BiexpResidualArray.mat','BiexpResidualArray','-v7.3');
save('YourPath\BiexpFitJacobianArray.mat','BiexpFitJacobianArray','-v7.3');
save('YourPath\BiexpFitOutputFlagArray.mat','BiexpFitOutputFlagArray','-v7.3');
save('YourPath\MonoexpFitArray.mat','MonoexpFitArray','-v7.3');
save('YourPath\MonoexpResidualArray.mat','MonoexpResidualArray','-v7.3');

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