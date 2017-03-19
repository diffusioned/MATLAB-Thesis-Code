% Chap3DisplayKurtFitOnBiexpSER.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap3DisplayKurtFitOnBiexpSER()

% Load noise free parameter values for data organization
load('U:\Code\MATLAB\Thesis Chapter\Biexponential\Monte Carlo Data\Thesis 50K NF Params\NoiseFreeParameters.mat');

% Load full fit kurtosis RSS
load('C:\Users\ccha4217\Documents\Aux Thesis Fits\Kurtosis Chapter\SNR 25 Run Kurtosis on Biexp Truth\KurtRSSArray.mat');
KurtRSSFFArray = KurtRSSArray;
clear KurtRSSArray;

% % Load 2SD kurtosis RSS
% load('C:\Users\ccha4217\Documents\Aux Thesis Fits\Kurtosis Chapter\SNR 25 Kurt Fit to Biexp 2SD\Kurt2SDRSSArray.mat');
% KurtRSS2SDArray = KurtRSSArray;
% clear KurtRSSArray;

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
MeanSERArray = zeros(SF1BinDim, D1D2BinDim);
MedianSERArray = zeros(SF1BinDim, D1D2BinDim);
StdDevSERArray = zeros(SF1BinDim, D1D2BinDim);
% Explicitly set the b-values and number of parameters here
BDim = 11;  TestParamDim = 3;
DOF = BDim - TestParamDim;

for i = 1:SF1BinDim
    for j = 1:D1D2BinDim
        CurBinIdxs = find(SF1Idxs == i & D1D2Idxs == j);
        % get SER for full fits
        CurRSSFF = KurtRSSFFArray(CurBinIdxs,:);
        SERFFArray = sqrt(CurRSSFF./DOF); % standard error of regression calculation
%         % and for 2SD fits
%         CurRSS2SD = KurtRSS2SDArray(CurBinIdxs,:);
%         SER2SDArray = sqrt(CurRSS2SD./DOF); % standard error of regression calculation
        
        MeanSERArray(i,j) = mean(SERFFArray(:));
        MedianSERArray(i,j) = median(SERFFArray(:));
        StdDevSERArray(i,j) = std(SERFFArray(:));
    end
    i
end

% Display difference in the two methods
PlotFontSize = 18;
% Colormap for Mean SNR
figure('Position', [0,0,800,600]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(MeanSERArray(:,:,1))');set(h,'edgecolor','none');
% caxis([-0.0025 0.005]);
% title('Mean SNR');
xlabel('True \itSF_{1}\rm\bf','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
% colormap(jet); % gives a good contrast for negative numbers
HandleC =  colorbar('location', 'EastOutside','Ticks',[0.04,0.05,0.06,0.07,0.08,0.09]);
HandleCLabel = ylabel(HandleC,'SER');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];
