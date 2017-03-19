% Chap4DisplaySelectionPctVsSkewedParameters.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap4DisplaySelectionPctVsSkewedParameters()

% Load noise free parameter values for data organization
load('YourPath\BiexpNoiseFreeParameters.mat');

% Load model fits
load('YourPath\BiexpFitArray11BValues.mat');
load('YourPath\MonoexpFitArray11BValues.mat');
load('YourPath\KurtFitArray11BValues.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load diagnostic metrics for biexp fits
load('YourPath\BiexpRSS11BVArray.mat');

% Kurtosis fits
load('YourPath\KurtRSS11BVArray.mat');

% Monoexp fits
load('YourPath\MonoexpRSS11BVArray.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare the AIC scores for biexp truth and display the plots
% Grab the dimensions
[NoiseFreeDim ~] = size(BiexpNFTestParameterArray);
[~, NoisySigDim, ~] = size(KurtRSS11BVArray);

% Explicitly set
BDim = 11;

% Sort the noise free parameters for display 
NFSF1ValueArray = BiexpNFTestParameterArray(:,1)./(BiexpNFTestParameterArray(:,1)+BiexpNFTestParameterArray(:,2));
NFD1D2RatioArray = BiexpNFTestParameterArray(:,3)./BiexpNFTestParameterArray(:,4);
SF1Dim = sqrt(length(NFSF1ValueArray));  D1D2Dim = sqrt(length(NFD1D2RatioArray)); % Only works if square array otherwise explicitly set

% Create AIC values, as well as AIC greater than arrays
MonoexpAICArray = zeros(NoiseFreeDim, NoisySigDim);
KurtAICArray = zeros(NoiseFreeDim, NoisySigDim);
BiexpAICArray = zeros(NoiseFreeDim, NoisySigDim);
BiexpAICMinusKurtosisAICArray = zeros(NoiseFreeDim, NoisySigDim);
BiexpAICMinusMonoexpAICArray = zeros(NoiseFreeDim, NoisySigDim);
KurtosisAICMinusMonoexpAICArray = zeros(NoiseFreeDim, NoisySigDim);

% Calculate the AIC values
MonoexpParam = 2+1; KurtParam = 3+1;  BiexpParam = 4+1; % The +1 is to include the variance needed with RSS fits
MonoexpAICArray = BDim.*log(MonoexpRSS11BVArray/BDim) + 2*(MonoexpParam);
KurtAICArray = BDim.*log(KurtRSS11BVArray/BDim) + 2*(KurtParam);
BiexpAICArray = BDim.*log(BiexpRSS11BVArray/BDim) + 2*(BiexpParam);
BiexpAICMinusKurtosisAICArray = BiexpAICArray - KurtAICArray;
BiexpAICMinusMonoexpAICArray = BiexpAICArray - MonoexpAICArray;
KurtosisAICMinusMonoexpAICArray = KurtAICArray - MonoexpAICArray;

% Create arrays for statistical calculations
MeanBiexpAICMinusKurtosisAICArray = zeros(SF1Dim, D1D2Dim);
MeanBiexpAICMinusMonoexpAICArray = zeros(SF1Dim, D1D2Dim);
MeanKurtosisAICMinusMonoexpAICArray = zeros(SF1Dim, D1D2Dim);
LowestAICScorePctArray  = zeros(SF1Dim, D1D2Dim, 3); % For 3 models

% Head to head models for selection illustration
BiexpAICSelectedVersusMonoexpArray = zeros(SF1Dim, D1D2Dim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Skewness Check
SkewnessCheckBiexpArray = zeros(SF1Dim,D1D2Dim,4);
SkewnessCheckKurtArray = zeros(SF1Dim,D1D2Dim,3);
SkewnessCheckMonoexpArray = zeros(SF1Dim,D1D2Dim,2);

% Lilliefors Test Check
LTestCheckBiexpArray = zeros(SF1Dim,D1D2Dim,4);
LTestCheckKurtArray = zeros(SF1Dim,D1D2Dim,3);
LTestCheckMonoexpArray = zeros(SF1Dim,D1D2Dim,2);
LTestCombinedBiexpArray = zeros(SF1Dim,D1D2Dim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop and calculate statisticals for each noise-free signal
for i = 1:length(NFSF1ValueArray)
    [sf,d] = ind2sub([SF1Dim, D1D2Dim],i);

    % Mean calculations
    CurBiexpAICMinusKurtosisAIC = BiexpAICMinusKurtosisAICArray(i,:);
    CurBiexpAICMinusMonoexpAIC = BiexpAICMinusMonoexpAICArray(i,:);
    CurKurtosisAICMinusMonoexpAIC = KurtosisAICMinusMonoexpAICArray(i,:);
    MeanBiexpAICMinusKurtosisAICArray(sf,d) = mean(CurBiexpAICMinusKurtosisAIC(:));
    MeanBiexpAICMinusMonoexpAICArray(sf,d) = mean(CurBiexpAICMinusMonoexpAIC(:));
    MeanKurtosisAICMinusMonoexpAICArray(sf,d) = mean(CurKurtosisAICMinusMonoexpAIC(:));

    % Find the lowest value
    AICScores = [MonoexpAICArray(i,:); KurtAICArray(i,:); BiexpAICArray(i,:);];
    [~,LowestAICScores] = min(AICScores);
    % Monoexp = 1, Kurtosis = 2, Biexp = 3;
    LowestAICScorePctArray(sf,d,1) = sum(LowestAICScores == 1)/NoisySigDim * 100.; % Monoexp
    LowestAICScorePctArray(sf,d,2) = sum(LowestAICScores == 2)/NoisySigDim * 100.; % Kurtosis
    LowestAICScorePctArray(sf,d,3) = sum(LowestAICScores == 3)/NoisySigDim * 100.; % Biexp
    
    % Head to head
    HeadToHeadAICBiexpVsMonoexp = [MonoexpAICArray(i,:); BiexpAICArray(i,:);];
    [~,LowestBiexpVsMonoexpAIC] = min(HeadToHeadAICBiexpVsMonoexp);
    BiexpAICSelectedVersusMonoexpArray(sf,d) = sum(LowestBiexpVsMonoexpAIC == 2)/NoisySigDim * 100;

    SkewnessCheckBiexpArray(sf,d,:) = (mean(BiexpFitArray(i,:,:))-median(BiexpFitArray(i,:,:)));
    SkewnessCheckKurtArray(sf,d,:) = (mean(KurtFitArray(i,:,:))-median(KurtFitArray(i,:,:)));
    SkewnessCheckMonoexpArray(sf,d,:) = (mean(MonoexpFitArray(i,:,:))-median(MonoexpFitArray(i,:,:)));

    % Lilliefors Test check
    AlphaSig = 0.001;
    LTestCheckBiexpArray(sf,d,1) = lillietest(squeeze(BiexpFitArray(i,:,1)),'Alpha',AlphaSig);
    LTestCheckBiexpArray(sf,d,2) = lillietest(squeeze(BiexpFitArray(i,:,2)),'Alpha',AlphaSig);
    LTestCheckBiexpArray(sf,d,3) = lillietest(squeeze(BiexpFitArray(i,:,3)),'Alpha',AlphaSig);
    LTestCheckBiexpArray(sf,d,4) = lillietest(squeeze(BiexpFitArray(i,:,4)),'Alpha',AlphaSig);
    LTestCheckKurtArray(sf,d,1) = lillietest(squeeze(KurtFitArray(i,:,1)),'Alpha',AlphaSig);
    LTestCheckKurtArray(sf,d,2) = lillietest(squeeze(KurtFitArray(i,:,2)),'Alpha',AlphaSig);
    LTestCheckKurtArray(sf,d,3) = lillietest(squeeze(KurtFitArray(i,:,3)),'Alpha',AlphaSig);
    LTestCheckMonoexpArray(sf,d,1) = lillietest(squeeze(MonoexpFitArray(i,:,1)),'Alpha',AlphaSig);
    LTestCheckMonoexpArray(sf,d,2) = lillietest(squeeze(MonoexpFitArray(i,:,2)),'Alpha',AlphaSig);

    % Check the combined values
    if LTestCheckBiexpArray(sf,d,1) == 0 || LTestCheckBiexpArray(sf,d,2) == 0 || LTestCheckBiexpArray(sf,d,3) == 0 || LTestCheckBiexpArray(sf,d,4) == 0
        LTestCombinedBiexpArray(sf,d) = 0;
    else
        LTestCombinedBiexpArray(sf,d) = 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hard coded from the parameter creation code
SF1Values = 0:1/69:1;
D1D2RatioValues = 2:18/69:20;
PlotFontSize = 12;

% Plot skewness values
figure('Position', [0,0,1000,800]);
subplot(2,2,1);h=pcolor(SF1Values, D1D2RatioValues, squeeze(SkewnessCheckBiexpArray(:,:,1))');caxis([-0.2 0.2]);set(h,'edgecolor','none');title('Skewness in \itA_{1}\rm\bf Estimates','fontsize', PlotFontSize+2);hold on;
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
HandleC = colorbar('location', 'EastOutside'); HandleCLabel = ylabel(HandleC,'Mean-Median'); drawnow; % MUST HAVE
LabelPos = HandleCLabel.Position; HandleCLabel.Rotation = -90; HandleCLabel.Position = [LabelPos(1) + 0.4, LabelPos(2), LabelPos(3)];

subplot(2,2,2);h=pcolor(SF1Values, D1D2RatioValues, squeeze(SkewnessCheckBiexpArray(:,:,2))');caxis([-0.2 0.2]);set(h,'edgecolor','none');title('Skewness in \itA_{2}\rm\bf Estimates','fontsize', PlotFontSize+2);hold on;
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
HandleC = colorbar('location', 'EastOutside'); HandleCLabel = ylabel(HandleC,'Mean-Median'); drawnow; % MUST HAVE
LabelPos = HandleCLabel.Position; HandleCLabel.Rotation = -90; HandleCLabel.Position = [LabelPos(1) + 0.4, LabelPos(2), LabelPos(3)];

subplot(2,2,3);h=pcolor(SF1Values, D1D2RatioValues, squeeze(SkewnessCheckBiexpArray(:,:,3))');caxis([-0.2 0.2]);set(h,'edgecolor','none');title('Skewness in \itD_{1}\rm\bf Estimates','fontsize', PlotFontSize+2);hold on;
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
HandleC = colorbar('location', 'EastOutside'); HandleCLabel = ylabel(HandleC,'Mean-Median'); drawnow; % MUST HAVE
LabelPos = HandleCLabel.Position; HandleCLabel.Rotation = -90; HandleCLabel.Position = [LabelPos(1) + 0.4, LabelPos(2), LabelPos(3)];

subplot(2,2,4);h=pcolor(SF1Values, D1D2RatioValues, squeeze(SkewnessCheckBiexpArray(:,:,4))');caxis([-0.2 0.2]);set(h,'edgecolor','none');title('Skewness in \itD_{2}\rm\bf Estimates','fontsize', PlotFontSize+2);hold on;
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
HandleC = colorbar('location', 'EastOutside'); HandleCLabel = ylabel(HandleC,'Mean-Median'); drawnow; % MUST HAVE
LabelPos = HandleCLabel.Position; HandleCLabel.Rotation = -90; HandleCLabel.Position = [LabelPos(1) + 0.4, LabelPos(2), LabelPos(3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lilliefors Test Results
figure('Position', [0,0,1000,800]);
subplot(2,2,1);h=pcolor(SF1Values, D1D2RatioValues, squeeze(LTestCheckBiexpArray(:,:,1))');set(h,'edgecolor','none');title('\itA_{1}\rm\bf Estimates','fontsize', PlotFontSize+2);hold on;
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);

subplot(2,2,2);h=pcolor(SF1Values, D1D2RatioValues, squeeze(LTestCheckBiexpArray(:,:,2))');set(h,'edgecolor','none');title('\itA_{2}\rm\bf Estimates','fontsize', PlotFontSize+2);hold on;
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);

subplot(2,2,3);h=pcolor(SF1Values, D1D2RatioValues, squeeze(LTestCheckBiexpArray(:,:,3))');set(h,'edgecolor','none');title('\itD_{1}\rm\bf Estimates','fontsize', PlotFontSize+2);hold on;
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);

subplot(2,2,4);h=pcolor(SF1Values, D1D2RatioValues, squeeze(LTestCheckBiexpArray(:,:,4))');set(h,'edgecolor','none');title('\itD_{2}\rm\bf Estimates','fontsize', PlotFontSize+2);hold on;
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Combined Lillefors Test
figure('Position', [0,0,800,600]);
h=pcolor(SF1Values, D1D2RatioValues, LTestCombinedBiexpArray(:,:)');set(h,'edgecolor','none');hold on;
set(gca,'FontSize',18,'FontWeight','bold'); % Fix axis label font size
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', 18);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', 18);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selection Percentages for delta AIC with overlaid normality test contour
figure('Position', [0,0,800,600]);
h=pcolor(SF1Values, D1D2RatioValues, MeanBiexpAICMinusMonoexpAICArray(:,:)');colorbar;set(h,'edgecolor','none');hold on;
set(gca,'FontSize',18,'FontWeight','bold'); % Fix axis label font size
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', 18);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', 18);
HandleC = colorbar('location', 'EastOutside'); HandleCLabel = ylabel(HandleC,'\Delta AIC (Biexp-Monoexp)'); drawnow; % MUST HAVE
LabelPos = HandleCLabel.Position; HandleCLabel.Rotation = -90; HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];

flt_size = 3; avg_flt_size = 3; pad_size = flt_size;% window size
AvgFilter = fspecial('average', [avg_flt_size avg_flt_size]);
LTestAverageContour = padarray(LTestCombinedBiexpArray(:,:)', [pad_size pad_size], 'replicate'); %zero gradient pad
LTestContourImage = imcrop(filter2(AvgFilter, LTestAverageContour), [pad_size + 1, pad_size + 1, length(SF1Values) - 1, length(D1D2RatioValues) - 1]);
contour(SF1Values, D1D2RatioValues, LTestContourImage(:,:), [1,1], 'LineColor', [0 0 0], 'LineWidth', 2);