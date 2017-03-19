% Chap4NLLSandLOOCVFit3ModelsToBiexpTruth.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap4NLLSandLOOCVFit3ModelsToBiexpTruth()

load('YourPath\BiexpNoisySignals_SNR25_11BValues.mat');
[NFSigDim NoisySigDim BDim] = size(BiexpNoisySignalArray);

% Set BValueArray back to itself, as cluster parfor can't see it - MATLAB strangeness???
TmpBValues = BValueArray;
BValueArray = TmpBValues;

% Calculate the noise std. dev.
TrueSTDNoise = 0.005 % For SNR of 25 = 0.04, SNR 200 = 0.005

% Create the array to store the fitted parameters based on the three spatial dimensions plus the axes and parameters.
% Make these all 1D for the parfor array and then reshape after.
MonoexpFitArray = zeros(NFSigDim, NoisySigDim, 2); % Amp, ADC
MonoexpResidualArray = zeros(NFSigDim, NoisySigDim, BDim);
MonoexpLOOCVScoreArray = zeros(NFSigDim,NoisySigDim);

BiexpFitArray = zeros(NFSigDim, NoisySigDim,4); % A1,A2,D1,D2
BiexpResidualArray = zeros(NFSigDim, NoisySigDim, BDim);
BiexpFitJacobianArray = zeros(NFSigDim, NoisySigDim,BDim,4); % Must be 4x4 to accommodate parameter Jacobians
BiexpFitOutputFlagArray = zeros(NFSigDim, NoisySigDim,1); % Saves the algorithm return code
BiexpLOOCVScoreArray = zeros(NFSigDim,NoisySigDim);

KurtFitArray = zeros(NFSigDim, NoisySigDim,3); % Amp, ADC, K
KurtResidualArray = zeros(NFSigDim, NoisySigDim, BDim);
KurtFitJacobianArray = zeros(NFSigDim, NoisySigDim, BDim,3); % Must be 4x4 to accommodate parameter Jacobians
KurtFitOutputFlagArray = zeros(NFSigDim, NoisySigDim,1); % Saves the algorithm return code
KurtLOOCVScoreArray = zeros(NFSigDim,NoisySigDim);

FitOpt = optimset('Display', 'off', 'Algorithm', 'trust-region-reflective', 'Jacobian', 'on');

% Now loop through and calculate all values for each voxel.  Flatten the 2D signal array to 1D.
for i = 1:NFSigDim
    parfor j = 1:NoisySigDim
        % FITTING
        VoxelSignal = squeeze(BiexpNoisySignalArray(i,j,:))';

        % Remove bias using Gudbjarsson equation
        % For each signal add a random amount of small error to take into account the fact that the std of the noise would not be known
        STDNoise = TrueSTDNoise + randn(1)*TrueSTDNoise*0.1;
        CorrVoxelSignal = sqrt(abs(VoxelSignal.^2 - STDNoise.^2));
        % or not
        % CorrVoxelSignal = VoxelSignal;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % MONOEXPONENTIAL FIT
        StartGuess = [max(CorrVoxelSignal), 1.5];% 2 fits
        LB = [0., 0.];    UB = [2*max(CorrVoxelSignal), 4.];
        [RetParams, RSS, Residuals]  =  lsqcurvefit(@MonoexpDecayWithJac, StartGuess, BValueArray, CorrVoxelSignal, LB, UB, FitOpt);
        MonoexpFitArray(i,j,:) = [RetParams(1), RetParams(2)];
        MonoexpResidualArray(i,j,:) = Residuals;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % BIEXPONENTIAL FIT
        LB = [0, 0, 0., 0.]; %[0., 0., 0., 0.]; % with constraints on linear
        UB = [max(CorrVoxelSignal), max(CorrVoxelSignal), 4., 4.]; % [GuessA0, 0.004, GuessA0, 0.004]; % with constraints on linear
        StartGuess = [0.5, 0.5, 1., 1/6.];

        [RetParams, RSS, Residuals,XFlag,Output,RetLamb,RetJac]  =  lsqcurvefit(@BiexpDecayIndAmpsWithJac, StartGuess, BValueArray, CorrVoxelSignal, LB, UB, FitOpt);
        if RetParams(3) < RetParams(4)
            % D2 is actually the higher diffusion value
            BiexpFitArray(i,j,:) = [RetParams(2), RetParams(1), RetParams(4), RetParams(3)];
        else
            BiexpFitArray(i,j,:) = [RetParams(1), RetParams(2), RetParams(3), RetParams(4)];
        end
        BiexpResidualArray(i,j,:) = Residuals;
        BiexpFitOutputFlagArray(i,j) = XFlag;
        BiexpFitJacobianArray(i,j,:,:) = full(RetJac);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % KURTOSIS FIT
        % Construct the starting guess and bounding vectors
        StartGuess = [max(CorrVoxelSignal), 1., 0.5];% 2 fits
        LB = [0, 0, -1.]; %[0., 0., 0., 0.]; % AMP, ADC, K
        UB = [2*max(CorrVoxelSignal), 4., 2.]; % AMP, ADC, K
        [RetParams, RSS, Residuals,XFlag,Output,RetLamb,RetJac]  =  lsqcurvefit(@KurtDecayWithJac, StartGuess, BValueArray, CorrVoxelSignal, LB, UB, FitOpt);
        KurtFitArray(i,j,:) = [RetParams(1), RetParams(2), RetParams(3)];
        KurtResidualArray(i,j,:) = Residuals;
        KurtFitJacobianArray(i,j,:,:) = full(RetJac);
        KurtFitOutputFlagArray(i,j,:) = XFlag;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % LEAVE-ONE-OUT CROSS-VALIDATION FITTING SELECTION
        % Setup the needed values beforehand
        LOOMonoexpTotalMSE = 0;
        LOOBiexpTotalMSE = 0;
        LOOKurtTotalMSE = 0;
        
        % Now loop through the b-values to test the LOOCV
        for cv = 1:BDim
            % Remove the current b-value and it's corresponding signal from the test set
            LOOBValues = BValueArray(BValueArray ~= BValueArray(cv));
            LOOSignal = VoxelSignal(VoxelSignal ~= VoxelSignal(cv));
            % FITTING - Remove bias using Gudbjarsson equation
            LOOSignal = sqrt(abs(LOOSignal.^2 - STDNoise.^2));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % MONOEXPONENTIAL FIT
            StartGuess = [max(LOOSignal), 1.5];% 2 fits
            LB = [0., 0.];    UB = [2*max(LOOSignal), 4];
            [RetParams, RSS, Residuals] = lsqcurvefit(@MonoexpDecayWithJac, StartGuess, LOOBValues, LOOSignal, LB, UB, FitOpt); % Jacobian on

            % Calculate monoexp loo sum
            YEst = MonoexpDecay(RetParams,BValueArray(cv));
            LOOMonoexpTotalMSE = LOOMonoexpTotalMSE + (YEst - VoxelSignal(cv)).^2;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % BIEXPONENTIAL FIT
            GuessA1 = 0.5*max(LOOSignal);  GuessA2  = 0.5*max(LOOSignal);
            GuessADC1 = 1.;  GuessADC2 = GuessADC1/6.;
            Guesses = [GuessA1, GuessA2, GuessADC1, GuessADC2];% 4 fits
            LB = [0, 0, 0., 0.];  UB = [1, 1, 4., 4.]; 
            [RetParams, RSS, Residuals] = lsqcurvefit(@BiexpDecayIndAmpsWithJac, Guesses, LOOBValues, LOOSignal, LB, UB, FitOpt); % Jacobian on

            % Calculate biexp loo sum
            YEst = BiexpDecay(RetParams,BValueArray(cv));
            LOOBiexpTotalMSE = LOOBiexpTotalMSE + (YEst - VoxelSignal(cv)).^2;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % KURTOSIS FIT
            StartGuess = [max(LOOSignal), 1., 0.5];% 2 fits
            LB = [0., 0., -1.];    UB = [2*max(LOOSignal), 4., 2.];
            [RetParams, RSS, Residuals] = lsqcurvefit(@KurtDecayWithJac, StartGuess, LOOBValues, LOOSignal, LB, UB, FitOpt); % Jacobian on

            % Calculate kurtosis loo sum
            YEst = KurtDecay(RetParams,BValueArray(cv));
            LOOKurtTotalMSE = LOOKurtTotalMSE + (YEst - VoxelSignal(cv)).^2;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
       
        MonoexpLOOCVScoreArray(i,j) = LOOMonoexpTotalMSE;
        BiexpLOOCVScoreArray(i,j) = LOOBiexpTotalMSE;
        KurtLOOCVScoreArray(i,j) = LOOKurtTotalMSE;
    end
    i
    datestr(now)
end

% Save all data in separate values
% Large files must be saved with the '-v7.3' flag
save('YourPath\MonoexpFitArray11BValues.mat','MonoexpFitArray','-v7.3');
save('YourPath\MonoexpResidualArray11BValues.mat','MonoexpResidualArray','-v7.3');
save('YourPath\MonoexpLOOCVScoreArray11BValues.mat','MonoexpLOOCVScoreArray','-v7.3');

save('YourPath\BiexpFitArray11BValues.mat','BiexpFitArray','-v7.3');
save('YourPath\BiexpResidualArray11BValues.mat','BiexpResidualArray','-v7.3');
save('YourPath\BiexpFitJacobianArray11BValues.mat','BiexpFitJacobianArray','-v7.3');
save('YourPath\BiexpFitOutputFlagArray11BValues.mat','BiexpFitOutputFlagArray','-v7.3');
save('YourPath\BiexpLOOCVScoreArray11BValues.mat','BiexpLOOCVScoreArray','-v7.3');

save('YourPath\KurtFitArray11BValues.mat','KurtFitArray','-v7.3');
save('YourPath\KurtResidualArray11BValues.mat','KurtResidualArray','-v7.3');
save('YourPath\KurtFitJacobianArray11BValues.mat','KurtFitJacobianArray','-v7.3');
save('YourPath\KurtFitOutputFlagArray11BValues.mat','KurtFitOutputFlagArray','-v7.3');
save('YourPath\KurtLOOCVScoreArray11BValues.mat','KurtLOOCVScoreArray','-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function signal = MonoexpDecay(params, xdata)
amp  = params(1); adc = params(2); 
signal = amp*exp(-adc*xdata);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function signal = BiexpDecay(params, xdata)
a1  = params(1); a2 = params(2); adc1 = params(3); adc2 = params(4);
signal = a1*exp(-adc1*xdata) + a2*exp(-adc2*xdata);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function signal = KurtDecay(params, xdata)
amp = params(1); adc = params(2); k = params(3);
signal = amp.*exp(-adc.*xdata + (k.*adc.^2.*xdata.^2)/6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signal,jacobian] = MonoexpDecayWithJac(params, xdata)
amp  = params(1); adc = params(2);
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
function [signal,jacobian] = KurtDecayWithJac(params, xdata)
amp = params(1); adc = params(2); k = params(3);
signal = amp.*exp(-adc.*xdata + (k.*adc.^2.*xdata.^2)/6);
% jacobian = (1/3.*amp.*adc).*(adc*k.*xdata-3).*exp(1/6.*adc.*xdata.*(adc*k.*xdata-6)); % This is d/dx
JacobPt1 = exp(1/6.*adc.^2*k.*xdata.^2-adc.*xdata); % This is d/d(amp)
JacobPt2 = (1/3.*amp.*xdata).*(adc*k.*xdata-3).*exp(1/6.*adc.*xdata.*(adc*k.*xdata-6)); % This is d/d(adc)
JacobPt3 = (1/6.*amp.*adc.^2.*xdata.^2).*exp(1/6.*adc.*xdata.*(adc*k.*xdata-6)); % This is d/k
jacobian = [JacobPt1' JacobPt2' JacobPt3']; % Complete Jacobian       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%