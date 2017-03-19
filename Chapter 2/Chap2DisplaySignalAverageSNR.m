% Chap2DisplaySignalAverageSNR.m
% MATLAB file for fitting simulated data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap2DisplaySignalAverageSNR()

% This is the main NLLS biexp fit function using the standard MATLAB lsqcurvefit function to fit the Monte Carlo simulated data for my thesis.
load('YourPath\NoiseFreeParameters.mat');
load('YourPath\NoiseFreeSignals.mat');

[NFSigDim BDim] = size(NFSignalArray);

% Calculate the noise std. dev.
TrueSTDNoise = 0.04 % For SNR of 25 = 0.04

% Create the array to store the fitted parameters based on the three spatial dimensions plus the axes and parameters.
% Make these all 1D for the parfor array and then reshape after.
MeanSNRArray = zeros(1, NFSigDim); % A1,A2,D1,D2,RSS
% For the max signal value to cut off the SNR
MeanSNR2SDArray = zeros(1, NFSigDim);

% Now loop through and calculate all values for each voxel.  Flatten the 2D signal array to 1D.
for i = 1:NFSigDim
    % Divide the noise-free measured signal at each b-value by the STD of the noise to get the SNR at each point
    SNRPointArray = NFSignalArray(i,:)/TrueSTDNoise;
    MeanSNRArray(i) = mean(SNRPointArray);
%     VoxelSignal = squeeze(NoisySignalArray(i,j,:))';
    MaxSignal = (NFSignalArray(i,:) > 2*TrueSTDNoise);
    % Find the first signal equal to zero and go back one
    [~,MinZeroIndex] = find(MaxSignal == 0);
    if(~isempty(MinZeroIndex))
        % Set to where the last signal was above 2SD
        MaxSignalIndex = min(MinZeroIndex)-1;
    end
    
    % Now limit the SNR point array to get the max
    MeanSNR2SDArray(i) = mean(SNRPointArray(1:MaxSignalIndex));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort the noise free parameters into bins - create 100x100 bins
SF1BinEdges = 0:0.01:1; SF1BinDim = 100;
D1D2BinEdges = 2:0.18:20; D1D2BinDim = 100;
EdgeVals = {SF1BinEdges, D1D2BinEdges}; % Cell array
NFSF1ValueArray = NFTestParameterArray(:,1)./(NFTestParameterArray(:,1)+NFTestParameterArray(:,2));
NFD1D2RatioArray = NFTestParameterArray(:,3)./NFTestParameterArray(:,4);
[SF1Cnt,~,SF1Idxs] = histcounts(NFSF1ValueArray,SF1BinEdges);
[D1D2Cnt,~,D1D2Idxs] = histcounts(NFD1D2RatioArray,D1D2BinEdges);
BinMeanSNRArray = zeros(SF1BinDim, D1D2BinDim);
Bin2SDMeanSNRArray = zeros(SF1BinDim, D1D2BinDim);

% Combine all errors into distribution metrics by parameter
for i = 1:SF1BinDim
    for j = 1:D1D2BinDim
        CurBinIdxs = find(SF1Idxs == i & D1D2Idxs == j);
        MeanSNRValues = MeanSNRArray(CurBinIdxs);
        MeanSNR2SDValues = MeanSNR2SDArray(CurBinIdxs);
        BinMeanSNRArray(i,j) = mean(MeanSNRValues);
        Bin2SDMeanSNRArray(i,j) = mean(MeanSNR2SDValues);


    end
end

PlotFontSize = 18;
% Colormap for Mean SNR
figure('Position', [0,0,800,600]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(BinMeanSNRArray(:,:,1))');set(h,'edgecolor','none');
% title('Mean SNR');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
% caxis([12 21]);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
% colormap(jet); % gives a good contrast for negative numbers
HandleC =  colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'SNR');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];

% or for the 2SD limited SNR
figure('Position', [0,0,800,600]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(Bin2SDMeanSNRArray(:,:,1))');set(h,'edgecolor','none');
% title('Mean SNR');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
% caxis([0 0.01]);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
% colormap(jet); % gives a good contrast for negative numbers
HandleC =  colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'SNR');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];
