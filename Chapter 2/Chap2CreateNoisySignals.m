% Chap2CreateNoisySignals.m
% MATLAB file for creating simulated data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap2CreateNoisySignals()
filename = 'YourPath\YourNoiseFreeSignals.mat'; mFile = matfile(filename);
NFSignalArray = mFile.NFSignalArray;
[NumberNFSignals ParamDim] = size(NFSignalArray);
BValueArray = mFile.BValueArray;

% This number sets the number of noisy signals you want to create for each noise-free signal
NumberNoisySignalsPerNFSignal = 200;
SNR = 25.;   STDNoise = 1/SNR;
BDim = length(BValueArray);
NoisySignalArray = zeros(NumberNFSignals, NumberNoisySignalsPerNFSignal, BDim);

% First loop through and create all noise free signals and for each NF signal add noise the specified number of times
for nf = 1:NumberNFSignals
    for s = 1:NumberNoisySignalsPerNFSignal
        for b = 1:length(BValueArray)
            NoisySignalArray(nf,s,b) = AddMagnitudeNoise(NFSignalArray(nf,b), STDNoise);
        end
    end
end

% Save in your desired directory with your desired name.  I ran this multiple times to create different sets
% with different SNR values, hence I add the SNR to the file name
save(sprintf('YourPath\YourNoisySignals_SNR%d.mat', floor(SNR)),'NoisySignalArray','BValueArray');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = AddMagnitudeNoise(a, s)
m = sqrt((a + s*randn(1)).^2 + (s*randn(1)).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%