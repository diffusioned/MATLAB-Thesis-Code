% Chap3CreateMonoexpNoiseFreeParameters.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap3CreateMonoexpNoiseFreeParameters()

NumberOfValues = 1000;
AmpValue = 1.;
MinADC = 0.05;   MaxADC = 1.;

NFTestParameterArray = zeros(NumberOfValues, 2); % Array parameters are A1, A2, D1, D2
ADCValues = (MaxADC-MinADC).*rand(NumberOfValues, 1) + MinADC;

NFTestParameterArray(:,1) = AmpValue; % Amplitude here is always 1
NFTestParameterArray(:,2) = ADCValues;

% Save the data
save('YourPath\MonoexpNFParametersForKurtosis.mat','NFTestParameterArray');