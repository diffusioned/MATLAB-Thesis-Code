% Chap3DisplayKurtosisFitJacobian.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap3DisplayKurtosisFitJacobian()

% Change the directories for displaying different data sets
% Load noise free parameter values for data organization
load('YourPath\NoiseFreeParameters.mat');

% Load kurt fits on biexp
load('YourPath\KurtFitArray.mat');

% Grab the dimensions
[NoiseFreeDim, ~] = size(NFTestParameterArray);
[~, NoisySigDim, ~] = size(KurtFitArray);
% Override the value of TestParamDim
TestParamDim = 3;

% Sort the noise free parameters into bins - create 100x100 bins
SF1BinEdges = 0:0.01:1; SF1BinDim = 100;
D1D2BinEdges = 2:0.18:20; D1D2BinDim = 100;
EdgeVals = {SF1BinEdges, D1D2BinEdges}; % Cell array
NFSF1ValueArray = NFTestParameterArray(:,1)./(NFTestParameterArray(:,1)+NFTestParameterArray(:,2));
NFD1D2RatioArray = NFTestParameterArray(:,3)./NFTestParameterArray(:,4);
[SF1Cnt,~,SF1Idxs] = histcounts(NFSF1ValueArray,SF1BinEdges);
[D1D2Cnt,~,D1D2Idxs] = histcounts(NFD1D2RatioArray,D1D2BinEdges);

PlotFontSize = 18;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Jacobian condition
load('YourPath\KurtJacCondArray.mat');
% Check some of the metrics
LargeCount = zeros(SF1BinDim, D1D2BinDim);
PctCount = zeros(SF1BinDim, D1D2BinDim);
% BinnedCondArray = 
CondNumThreshold = 20;
for i = 1:SF1BinDim
    for j = 1:D1D2BinDim
        CurBinIdxs = find(SF1Idxs == i & D1D2Idxs == j);
        CurCondArray = KurtJacCondArray(CurBinIdxs,:);
        LargeCount(i,j) = length(find(CurCondArray(:) > CondNumThreshold));
        PctCount(i,j) = LargeCount(i,j)/length(CurBinIdxs);
    end
end

figure('Position', [0,0,1000,800]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(PctCount(:,:))');caxis([0 100]);
set(h,'edgecolor','none');
% title('Number of regression fits with Jacobian condition greater than 100');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% caxis([0.014 0.028]); % SNR 25
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'EastOutside');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%