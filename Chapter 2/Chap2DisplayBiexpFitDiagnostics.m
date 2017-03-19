% Chap2DisplayBiexpFitDiagnostics.m
% MATLAB file for display fit diagnostics in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap2DisplayBiexpFitDiagnostics()

% Change the directories for displaying different data sets
% Load noise free parameter values for data organization
load('YourPath\YourNoiseFreeParameters.mat');

% Load biexp fits
load('YourPath\BiexpFitArray.mat');

% Grab the dimensions
[NoiseFreeDim TestParamDim] = size(NFTestParameterArray);
[~, NoisySigDim, ~] = size(BiexpFitArray);

% Sort the noise free parameters into bins - create 100x100 bins
SF1BinEdges = 0:0.01:1; SF1BinDim = 100;
D1D2BinEdges = 2:0.18:20; D1D2BinDim = 100;
EdgeVals = {SF1BinEdges, D1D2BinEdges}; % Cell array
NFSF1ValueArray = NFTestParameterArray(:,1)./(NFTestParameterArray(:,1)+NFTestParameterArray(:,2));
NFD1D2RatioArray = NFTestParameterArray(:,3)./NFTestParameterArray(:,4);
[SF1Cnt,~,SF1Idxs] = histcounts(NFSF1ValueArray,SF1BinEdges);
[D1D2Cnt,~,D1D2Idxs] = histcounts(NFD1D2RatioArray,D1D2BinEdges);

PlotFontSize = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Jacobian condition
load('YourPath\JacCondArray.mat');
% [NFDim SigDim] = size(JacCondArray);
% JacobianCondArray = zeros(NFDim, SigDim,ParamDim,ParamDim);
% Check some of the metrics
MeanCondArray = zeros(SF1BinDim, D1D2BinDim);
StdDevArray = zeros(SF1BinDim, D1D2BinDim);
MedianCondArray = zeros(SF1BinDim, D1D2BinDim);
IQRCondArray = zeros(SF1BinDim, D1D2BinDim);
LargeCount = zeros(SF1BinDim, D1D2BinDim);
PctCount = zeros(SF1BinDim, D1D2BinDim);
% BinnedCondArray = 
CondNumThreshold = 1e6; % CHANGE THIS THRESHOLD TO PLOT A DIFFERENT VALUE
for i = 1:SF1BinDim
    for j = 1:D1D2BinDim
        CurBinIdxs = find(SF1Idxs == i & D1D2Idxs == j);
        CurCondArray = JacCondArray(CurBinIdxs,:);
        MeanCondArray(i,j) = mean(CurCondArray(:));
        StdDevArray(i,j) = std(CurCondArray(:));
        MedianCondArray(i,j) = median(CurCondArray(:));
        IQRCondArray(i,j) = iqr(CurCondArray(:));
        LargeCount(i,j) = length(find(CurCondArray(:) > CondNumThreshold));
        PctCount(i,j) = 100*LargeCount(i,j)/length(CurCondArray(:));
    end
end
% Group the noise free parameters into bins and translate the values
PlotFontSize = 18;
figure('Position', [0,0,800,600]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(PctCount(:,:))');
caxis([0 30]);set(h,'edgecolor','none');
% title('Mean Std Err of Regression');
xlabel('True \itSF_{1}\rm\bf','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC = colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'%');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter standard errors (from covariance matrix) and estimated confidence intervals
load('YourPath\CovMxArray.mat');
load('YourPath\BiexpFitArray.mat');

% Check some of the metrics
MeanParamSDArray = zeros(SF1BinDim, D1D2BinDim, TestParamDim);
MedianParamSDArray = zeros(SF1BinDim, D1D2BinDim, TestParamDim);
CIReliabilityArray = zeros(SF1BinDim, D1D2BinDim, TestParamDim); % Array to count times true value falls in CI
DOFModel = 11-4; % signals - fit parameters
alpha = 0.05; % 95% CI
TFactor = tinv(1-alpha/2,DOFModel);

for i = 1:SF1BinDim
    for j = 1:D1D2BinDim
        CurBinIdxs = find(SF1Idxs == i & D1D2Idxs == j);
        CurSEAmp1Array = sqrt(CovMxArray(CurBinIdxs,:,1,1));
        CurSEAmp2Array = sqrt(CovMxArray(CurBinIdxs,:,2,2));
        CurSEADC1Array = sqrt(CovMxArray(CurBinIdxs,:,3,3));
        CurSEADC2Array = sqrt(CovMxArray(CurBinIdxs,:,4,4));
        
        MeanParamSDArray(i,j,1) = mean(CurSEAmp1Array(:));
        MeanParamSDArray(i,j,2) = mean(CurSEAmp2Array(:));
        MeanParamSDArray(i,j,3) = mean(CurSEADC1Array(:));
        MeanParamSDArray(i,j,4) = mean(CurSEADC2Array(:));
        
        MedianParamSDArray(i,j,1) = median(CurSEAmp1Array(:));
        MedianParamSDArray(i,j,2) = median(CurSEAmp2Array(:));
        MedianParamSDArray(i,j,3) = median(CurSEADC1Array(:));
        MedianParamSDArray(i,j,4) = median(CurSEADC2Array(:));
        
        % Combine with the fit values for CI
        CurMinAmp1CIArray = BiexpFitArray(CurBinIdxs,:,1) - CurSEAmp1Array*TFactor;
        CurMinAmp2CIArray = BiexpFitArray(CurBinIdxs,:,2) - CurSEAmp2Array*TFactor;
        CurMinADC1CIArray = BiexpFitArray(CurBinIdxs,:,3) - CurSEADC1Array*TFactor;
        CurMinADC2CIArray = BiexpFitArray(CurBinIdxs,:,4) - CurSEADC2Array*TFactor;
        
        CurMaxAmp1CIArray = BiexpFitArray(CurBinIdxs,:,1) + CurSEAmp1Array*TFactor;
        CurMaxAmp2CIArray = BiexpFitArray(CurBinIdxs,:,2) + CurSEAmp2Array*TFactor;
        CurMaxADC1CIArray = BiexpFitArray(CurBinIdxs,:,3) + CurSEADC1Array*TFactor;
        CurMaxADC2CIArray = BiexpFitArray(CurBinIdxs,:,4) + CurSEADC2Array*TFactor;
        
        % Now count how many of the confidence intervals encompass the true value
%         NumSignalValues = length(CurBinIdxs)*NoisySigDim;
        InCICount = zeros(TestParamDim,1);
        for nf = 1:length(CurBinIdxs)
            CurNFIdx = CurBinIdxs(nf);
            for s = 1:NoisySigDim
                if(NFTestParameterArray(CurNFIdx,1) <= CurMaxAmp1CIArray(nf,s) &&  NFTestParameterArray(CurNFIdx,1) >= CurMinAmp1CIArray(nf,s))
                    InCICount(1) = InCICount(1) + 1;
                end
                if(NFTestParameterArray(CurNFIdx,2) <= CurMaxAmp2CIArray(nf,s) &&  NFTestParameterArray(CurNFIdx,2) >= CurMinAmp2CIArray(nf,s))
                    InCICount(2) = InCICount(2) + 1;
                end
                if(NFTestParameterArray(CurNFIdx,3) <= CurMaxADC1CIArray(nf,s) &&  NFTestParameterArray(CurNFIdx,3) >= CurMinADC1CIArray(nf,s))
                    InCICount(3) = InCICount(3) + 1;
                end               
                if(NFTestParameterArray(CurNFIdx,4) <= CurMaxADC2CIArray(nf,s) &&  NFTestParameterArray(CurNFIdx,4) >= CurMinADC2CIArray(nf,s))
                    InCICount(4) = InCICount(4) + 1;
                end

            end
        end
        
        TotalCount = length(CurBinIdxs)*NoisySigDim;
        CIReliabilityArray(i,j,:) = InCICount(:)/TotalCount * 100.;

    end
    i
end

PlotFontSize = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean Param SD
figure('Position', [0,0,1000,800]);
subplot(2,2,1);h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(MeanParamSDArray(:,:,1))');caxis([0 1000]);
set(h,'edgecolor','none');title('Mean \itA_{1}\rm\bf Parameter SD');
xlabel('True \itSF_{1}\rm\bf','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'SD','FontSize',PlotFontSize,'FontWeight','bold');
% Optional code to rotate the colorbar label 180 degrees
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 0.6, LabelPos(2), LabelPos(3)];
subplot(2,2,2);h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(MeanParamSDArray(:,:,2))');caxis([0 1000]);
set(h,'edgecolor','none');title('Mean \itA_{2}\rm\bf Parameter SD');
xlabel('True \itSF_{1}\rm\bf','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'SD','FontSize',PlotFontSize,'FontWeight','bold');
% Optional code to rotate the colorbar label 180 degrees
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 0.6, LabelPos(2), LabelPos(3)];
subplot(2,2,3);h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(MeanParamSDArray(:,:,3))');caxis([0 1000]);
set(h,'edgecolor','none');title('Mean \itD_{1}\rm\bf Parameter SD');
xlabel('True \itSF_{1}\rm\bf','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'SD','FontSize',PlotFontSize,'FontWeight','bold');
% Optional code to rotate the colorbar label 180 degrees
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 0.6, LabelPos(2), LabelPos(3)];
subplot(2,2,4);h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(MeanParamSDArray(:,:,4))');caxis([0 1000]);
set(h,'edgecolor','none');title('Mean \itD_{2}\rm\bf Parameter SD');
xlabel('True \itSF_{1}\rm\bf','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'SD','FontSize',PlotFontSize,'FontWeight','bold');
% Optional code to rotate the colorbar label 180 degrees
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 0.6, LabelPos(2), LabelPos(3)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Median SD
figure('Position', [0,0,1000,800]);
subplot(2,2,1);h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(MedianParamSDArray(:,:,1))');caxis([0 0.2]);
set(h,'edgecolor','none');title('Median \itA_{1}\rm\bf Parameter SD');
xlabel('True \itSF_{1}\rm\bf','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'SD','FontSize',PlotFontSize,'FontWeight','bold');
% Optional code to rotate the colorbar label 180 degrees
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 0.6, LabelPos(2), LabelPos(3)];
subplot(2,2,2);h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(MedianParamSDArray(:,:,2))');caxis([0 0.2]);
set(h,'edgecolor','none');title('Median \itA_{2}\rm\bf Parameter SD');
xlabel('True \itSF_{1}\rm\bf','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'SD','FontSize',PlotFontSize,'FontWeight','bold');
% Optional code to rotate the colorbar label 180 degrees
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 0.6, LabelPos(2), LabelPos(3)];
subplot(2,2,3);h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(MedianParamSDArray(:,:,3))');caxis([0 0.2]);
set(h,'edgecolor','none');title('Median \itD_{1}\rm\bf Parameter SD');
xlabel('True \itSF_{1}\rm\bf','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'SD','FontSize',PlotFontSize,'FontWeight','bold');
% Optional code to rotate the colorbar label 180 degrees
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 0.6, LabelPos(2), LabelPos(3)];
subplot(2,2,4);h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(MedianParamSDArray(:,:,4))');caxis([0 0.2]);
set(h,'edgecolor','none');title('Median \itD_{2}\rm\bf Parameter SD');
xlabel('True \itSF_{1}\rm\bf','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'SD','FontSize',PlotFontSize,'FontWeight','bold');
% Optional code to rotate the colorbar label 180 degrees
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 0.6, LabelPos(2), LabelPos(3)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 95% CI reliability values
figure('Position', [0,0,1000,1000]);
subplot(2,2,1);h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(CIReliabilityArray(:,:,1))');caxis([40 100]);
set(h,'edgecolor','none');title('Reliability of \itA_{1}\rm\bf 95% CI');
xlabel('True \itSF_{1}\rm\bf','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'% Within Interval','FontSize',PlotFontSize,'FontWeight','bold');
% Optional code to rotate the colorbar label 180 degrees
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 0.6, LabelPos(2), LabelPos(3)];
subplot(2,2,2);h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(CIReliabilityArray(:,:,2))');caxis([40 100]);
set(h,'edgecolor','none');title('Reliability of \itA_{2}\rm\bf 95% CI');
xlabel('True \itSF_{1}\rm\bf','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'% Within Interval','FontSize',PlotFontSize,'FontWeight','bold');
% Optional code to rotate the colorbar label 180 degrees
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 0.6, LabelPos(2), LabelPos(3)];
subplot(2,2,3);h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(CIReliabilityArray(:,:,3))');caxis([40 100]);
set(h,'edgecolor','none');title('Reliability of \itD_{1}\rm\bf 95% CI');
xlabel('True \itSF_{1}\rm\bf','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'% Within Interval','FontSize',PlotFontSize,'FontWeight','bold');
% Optional code to rotate the colorbar label 180 degrees
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 0.6, LabelPos(2), LabelPos(3)];
subplot(2,2,4);h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(CIReliabilityArray(:,:,4))');caxis([40 100]);
set(h,'edgecolor','none');title('Reliability of \itD_{2}\rm\bf 95% CI');
xlabel('True \itSF_{1}\rm\bf','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC =  colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'% Within Interval','FontSize',PlotFontSize,'FontWeight','bold');
% Optional code to rotate the colorbar label 180 degrees
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 0.6, LabelPos(2), LabelPos(3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and display the Biexponential Standard Error of Regression
load('YourPath\BiexpResidualArray.mat');

% Now group into bins and plot the mean and std dev of the standard error of regression
MeanSERArray = zeros(SF1BinDim, D1D2BinDim);
MedianSERArray = zeros(SF1BinDim, D1D2BinDim);
StdDevSERArray = zeros(SF1BinDim, D1D2BinDim);
KSNormalTestArray = zeros(SF1BinDim, D1D2BinDim);
% Extract the number of b-values (residuals)
[~,~,BDim] = size(BiexpResidualArray);
DOF = BDim - TestParamDim;

for i = 1:SF1BinDim
    for j = 1:D1D2BinDim
        CurBinIdxs = find(SF1Idxs == i & D1D2Idxs == j);
        % get all residuals for those fits
        CurSERArray = BiexpResidualArray(CurBinIdxs,:,:);
        CurRSS = sum(CurSERArray.^2,3);
        SERArray = sqrt(CurRSS./DOF); % standard error of regression calculation
        MeanSERArray(i,j) = mean(SERArray(:));
        MedianSERArray(i,j) = median(SERArray(:));
        StdDevSERArray(i,j) = std(SERArray(:));
        
        % NORMALITY TEST FOR RESIDUALS
        if not(isempty(CurSERArray(:)))
            KSNormalTestArray(i,j) = kstest(CurSERArray(:));
        end
    end
    i
end

% Thesis plot mean SER
PlotFontSize = 18;
figure('Position', [100,0,1000,800]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(MeanSERArray(:,:))');
caxis([0.036 0.042]);
colorbar('location', 'EastOutside');
set(h,'edgecolor','none');
% title('Mean Std Err of Regression');
xlabel('True \itSF_{1}\rm\bf','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

% Histogram of values
figure('Position', [0,0,1000,800]);
histogram(MeanSERArray(:),100);
xlabel('SER','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Bin Count','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find an "uncertain" signal and a "less uncertain" one
% Uncertain, Index: 111 (EXAMPLE), values: 0.0081    0.9919    1.0000    0.4850
% Less Uncertain, Index: 222 (EXAMPLE),  values: 0.5011    0.4989    1.0000    0.0504
UncertainResids = squeeze(BiexpResidualArray(111,:,:));
CertainResids = squeeze(BiexpResidualArray(222,:,:));
figure('Position', [100,100,1200,600]);
subplot(1,2,1);qqplot(squeeze(UncertainResids(:)));title('True Parameters \itSF_{1}\rm\bf = 0.01, \itD_{1} / D_{2}\rm\bf = 2.06');
set(gca,'FontSize',14, 'FontWeight', 'bold');
subplot(1,2,2);qqplot(squeeze(CertainResids(:)));title('True Parameters \itSF_{1}\rm\bf = 0.5, \itD_{1} / D_{2}\rm\bf = 19.85');
set(gca,'FontSize',14, 'FontWeight', 'bold');
% xlabel('log Signal Amplitude', 'fontsize', 24, 'fontweight', 'bold');
% ylabel('Diffusion Weighting (AU)', 'fontsize', 24, 'fontweight', 'bold');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
