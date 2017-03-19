% Chap3CompareFulland2SDFitOnKurtSER.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap3CompareFulland2SDFitOnKurtSER()

% Load noise free parameter values for data organization
load('YourPath\NoiseFreeParameters.mat');

% Load full fit kurtosis RSS
load('YourPath\KurtRSSArray.mat');
KurtRSSFFArray = KurtRSSArray;
clear KurtRSSArray;

% Load 2SD kurtosis RSS
load('YourPath\Kurt2SDRSSArray.mat');
KurtRSS2SDArray = KurtRSSArray;
clear KurtRSSArray;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sort the noise free parameters into bins - create 100x100 bins
SF1BinEdges = 0:0.01:1; SF1BinDim = 100;
D1D2BinEdges = 2:0.18:20; D1D2BinDim = 100;
EdgeVals = {SF1BinEdges, D1D2BinEdges}; % Cell array
NFSF1ValueArray = NFTestParameterArray(:,1)./(NFTestParameterArray(:,1)+NFTestParameterArray(:,2));
NFD1D2RatioArray = NFTestParameterArray(:,3)./NFTestParameterArray(:,4);
[SF1Cnt,~,SF1Idxs] = histcounts(NFSF1ValueArray,SF1BinEdges);
[D1D2Cnt,~,D1D2Idxs] = histcounts(NFD1D2RatioArray,D1D2BinEdges);

% Load and display the monoexponential Standard Error of Regression (SER)
% Now group into bins and plot the mean and std dev of the standard error of regression
MeanDiffSERArray = zeros(SF1BinDim, D1D2BinDim);
MedianDiffSERArray = zeros(SF1BinDim, D1D2BinDim);
StdDevDiffSERArray = zeros(SF1BinDim, D1D2BinDim);
% Explicitly set the b-values and number of parameters here
BDim = 11;  TestParamDim = 3;
DOF = BDim - TestParamDim;

for i = 1:SF1BinDim
    for j = 1:D1D2BinDim
        CurBinIdxs = find(SF1Idxs == i & D1D2Idxs == j);
        % get SER for full fits
        CurRSSFF = KurtRSSFFArray(CurBinIdxs,:);
        SERFFArray = sqrt(CurRSSFF./DOF); % standard error of regression calculation
        % and for 2SD fits
        CurRSS2SD = KurtRSS2SDArray(CurBinIdxs,:);
        SER2SDArray = sqrt(CurRSS2SD./DOF); % standard error of regression calculation
        
        MeanDiffSERArray(i,j) = mean(SER2SDArray(:)-SERFFArray(:));
        MedianDiffSERArray(i,j) = median(SER2SDArray(:)-SERFFArray(:));
        StdDevDiffSERArray(i,j) = std(SER2SDArray(:)-SERFFArray(:));
    end
end

% Display difference in the two methods
PlotFontSize = 18;
% Colormap for Mean SNR
figure('Position', [0,0,800,600]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(MeanDiffSERArray(:,:,1))');set(h,'edgecolor','none');
xlabel('True \itSF_{1}\rm\bf','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC = colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'Difference in SER');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];
