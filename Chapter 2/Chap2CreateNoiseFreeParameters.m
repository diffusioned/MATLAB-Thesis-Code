% Chap2LinePlotOfSignalRange.m
% MATLAB file for creating simulated data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap2CreateNoiseFreeParameters()
% This file creates random noise-free parameters for a biexponential signal test set
NumberOfValues = 50000; % How many samples you want
MinSF1Value = 0.;  MaxSF1Value = 1.; % The range for SF1
MinD1D2Ratio = 2.;   MaxD1D2Ratio = 20.; % The range for the ratio of D1/D2
NFTestParameterArray = zeros(NumberOfValues, 4); % Array parameters are A1, A2, D1, D2
SignalFraction1Values = (MaxSF1Value-MinSF1Value).*rand(NumberOfValues, 1) + MinSF1Value; % SF1 = A1 here because TrueD1Value = 1
Amp2Values = 1 - SignalFraction1Values;
D1D2RatioValues = (MaxD1D2Ratio-MinD1D2Ratio).*rand(NumberOfValues, 1) + MinD1D2Ratio;
D2Values = 1./D1D2RatioValues;

NFTestParameterArray(:,1) = SignalFraction1Values;
NFTestParameterArray(:,2) = Amp2Values;
NFTestParameterArray(:,3) = 1;
NFTestParameterArray(:,4) = D2Values;

% Save the data in a desired directory and call it what you want
save('YourPath\YourNoiseFreeParameters.mat','NFTestParameterArray');
