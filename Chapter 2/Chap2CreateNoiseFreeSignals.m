% Chap2LinePlotOfSignalRange.m
% MATLAB file for creating simulated data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap2CreateNoiseFreeSignals()
% Append your path to the filename below and change the .mat file if you renamed it something else.
filename = 'YourPath\YourNoiseFreeParameters.mat'; mFile = matfile(filename);
NFTestParameterArray = mFile.NFTestParameterArray;
[NumberNFSignals ParamDim] = size(NFTestParameterArray);

% Now we need the b-value array to be set here.  These normalized values from my thesis were chosen so that the amplitude of D1
% decreases to approximate 1/e times the original signal value.  There are also five values covering D1, one
% at the transition point where D1 is near zero and then five covering D2 mostly.
BValueArray = [0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6];
BDim = length(BValueArray);
NFSignalArray = zeros(NumberNFSignals, BDim);
% Now loop and create all signals
for nf = 1:NumberNFSignals
    NFSignalArray(nf,:) = CreateNoiseFreeSignal(NFTestParameterArray(nf,1), NFTestParameterArray(nf,2), NFTestParameterArray(nf,3), ...
                                                NFTestParameterArray(nf,4), BValueArray);
end
% Save the data in the directory you want with the file name you want
save('YourPath\YourNoiseFreeSignals.mat','NFSignalArray','BValueArray');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = CreateNoiseFreeSignal(a1, a2, adc1, adc2, b)
S = a1.*exp(-adc1.*b) + a2.*exp(-adc2.*b);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
