% Chap2CompareFulland2SDFitOnMonoSER.m
% MATLAB file for fitting simulated data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap2CompareFulland2SDFitOnMonoSER()

% Load noise free parameter values for data organization
load('YourPath\NoiseFreeParameters.mat');

% Load full fit monoexp RSS
load('YourPath\MonoexpRSSArray.mat');
MonoexpRSSFFArray = MonoexpRSSArray;
clear MonoexpRSSArray;

% Load monoexp RSS
load('YourPath\MonoexpRSSArray.mat');
MonoexpRSS2SDArray = MonoexpRSSArray;
clear MonoexpRSSArray;

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
BDim = 11;  TestParamDim = 2;
DOF = BDim - TestParamDim;

for i = 1:SF1BinDim
    for j = 1:D1D2BinDim
        CurBinIdxs = find(SF1Idxs == i & D1D2Idxs == j);
        % get SER for full fits
        CurRSSFF = MonoexpRSSFFArray(CurBinIdxs,:);
        SERFFArray = sqrt(CurRSSFF./DOF); % standard error of regression calculation
        % and for 2SD fits
        CurRSS2SD = MonoexpRSS2SDArray(CurBinIdxs,:);
        SER2SDArray = sqrt(CurRSS2SD./DOF); % standard error of regression calculation
        
        MeanDiffSERArray(i,j) = mean(SERFFArray(:)-SER2SDArray(:));
        MedianDiffSERArray(i,j) = median(SERFFArray(:)-SER2SDArray(:));
        StdDevDiffSERArray(i,j) = std(SERFFArray(:)-SER2SDArray(:));
    end
end

% Display difference in the two methods
PlotFontSize = 18;
% Colormap for Mean SNR
figure('Position', [0,0,800,600]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(MeanDiffSERArray(:,:,1))');set(h,'edgecolor','none');
caxis([0 0.01]);
% title('Mean SNR');
xlabel('True \itSF_{1}\rm\bf','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
% colormap(jet); % gives a good contrast for negative numbers
HandleC =  colorbar('location', 'EastOutside','Ticks',[0,0.002,0.004,0.006,0.008,0.01]);
HandleCLabel = ylabel(HandleC,'SER Difference');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1, LabelPos(2), LabelPos(3)];
