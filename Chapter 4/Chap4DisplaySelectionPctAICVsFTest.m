% Chap4DisplaySelectionPctAICVsFTest.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap4DisplaySelectionPctAICVsFTest()

% Load noise free parameter values for data organization
load('YourPath\BiexpNoiseFreeParameters.mat');

% Load model fits
load('YourPath\BiexpFitArray11BValues.mat');
load('YourPath\MonoexpFitArray11BValues.mat');
load('YourPath\KurtFitArray11BValues.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load diagnostic metrics for biexp fits
load('YourPath\BiexpRSS11BVArray.mat');
load('YourPath\BiexpLOOCVScoreArray11BValues.mat');

% Kurtosis fits
load('YourPath\KurtRSS11BVArray.mat');
load('YourPath\KurtLOOCVScoreArray11BValues.mat');

% Monoexp fits
load('YourPath\MonoexpRSS11BVArray.mat');
load('YourPath\MonoexpLOOCVScoreArray11BValues.mat');
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
BiexpAICcSelectedVersusMonoexpArray = zeros(SF1Dim, D1D2Dim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AICc - repeat for corrected AIC values
% Create AICc values, as well as AICc difference arrays
MonoexpAICcArray = zeros(NoiseFreeDim, NoisySigDim);
KurtAICcArray = zeros(NoiseFreeDim, NoisySigDim);
BiexpAICcArray = zeros(NoiseFreeDim, NoisySigDim);
BiexpAICcMinusKurtosisAICcArray = zeros(NoiseFreeDim, NoisySigDim);
BiexpAICcMinusMonoexpAICcArray = zeros(NoiseFreeDim, NoisySigDim);
KurtosisAICMinusMonoexpAICcArray = zeros(NoiseFreeDim, NoisySigDim);

% Calculate the AICc values
MonoexpAICcArray = MonoexpAICArray + 2*MonoexpParam*(MonoexpParam+1)/(BDim-MonoexpParam-1); % 11 B values = 3.4286
KurtAICcArray = KurtAICArray + 2*KurtParam*(KurtParam+1)/(BDim-KurtParam-1); % 11 B values = 6.6667
BiexpAICcArray = BiexpAICArray + 2*BiexpParam*(BiexpParam+1)/(BDim-BiexpParam-1); % 11 B values = 12
BiexpAICcMinusKurtosisAICcArray = BiexpAICcArray - KurtAICcArray;
BiexpAICcMinusMonoexpAICcArray = BiexpAICcArray - MonoexpAICcArray;
KurtosisAICcMinusMonoexpAICcArray = KurtAICcArray - MonoexpAICcArray;

% Create arrays for statistical calculations
MeanBiexpAICcMinusKurtosisAICcArray = zeros(SF1Dim, D1D2Dim);
MeanBiexpAICcMinusMonoexpAICcArray = zeros(SF1Dim, D1D2Dim);
MeanKurtosisAICcMinusMonoexpAICcArray = zeros(SF1Dim, D1D2Dim);
LowestAICcScorePctArray  = zeros(SF1Dim, D1D2Dim, 3); % For 3 models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assessment of LOOCV scores
BiexpLOOCVMinusKurtosisLOOCVArray = zeros(SF1Dim, D1D2Dim);
BiexpLOOCVMinusMonoexpLOOCVArray = zeros(SF1Dim, D1D2Dim);
KurtosisLOOCVMinusMonoexpLOOCVArray = zeros(SF1Dim, D1D2Dim);
BiexpLOOCVMinusKurtosisLOOCVArray = BiexpLOOCVScoreArray - KurtLOOCVScoreArray;
BiexpLOOCVMinusMonoexpLOOCVArray = BiexpLOOCVScoreArray - MonoexpLOOCVScoreArray;
KurtosisLOOCVMinusMonoexpLOOCVArray = KurtLOOCVScoreArray - MonoexpLOOCVScoreArray;
% Stats - NOTE LOOCV needs median because several values of the kurtosis LOOCV are near 10^17
MeanBiexpLOOCVMinusKurtosisLOOCVArray = zeros(SF1Dim, D1D2Dim);
MeanBiexpLOOCVMinusMonoexpLOOCVArray = zeros(SF1Dim, D1D2Dim);
MeanKurtosisLOOCVMinusMonoexpLOOCVArray = zeros(SF1Dim, D1D2Dim);
LowestLOOCVScorePctArray  = zeros(SF1Dim, D1D2Dim, 3); % For 3 models

% FTestScoreArray = zeros(SF1Dim, D1D2Dim);
FTestBiexpBetterPctArray = zeros(SF1Dim, D1D2Dim);

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
    % AICc
    CurBiexpAICcMinusKurtosisAICc = BiexpAICcMinusKurtosisAICcArray(i,:);
    CurBiexpAICcMinusMonoexpAICc = BiexpAICcMinusMonoexpAICcArray(i,:);
    CurKurtosisAICcMinusMonoexpAICc = KurtosisAICcMinusMonoexpAICcArray(i,:);
    MeanBiexpAICcMinusKurtosisAICcArray(sf,d) = mean(CurBiexpAICcMinusKurtosisAICc(:));
    MeanBiexpAICcMinusMonoexpAICcArray(sf,d) = mean(CurBiexpAICcMinusMonoexpAICc(:));
    MeanKurtosisAICcMinusMonoexpAICcArray(sf,d) = mean(CurKurtosisAICcMinusMonoexpAICc(:));
    % LOOCV
    CurBiexpLOOCVMinusKurtosisLOOCV = BiexpLOOCVMinusKurtosisLOOCVArray(i,:);
    CurBiexpLOOCVMinusMonoexpLOOCV = BiexpLOOCVMinusMonoexpLOOCVArray(i,:);
    CurKurtosisLOOCVMinusMonoexpLOOCV = KurtosisLOOCVMinusMonoexpLOOCVArray(i,:);
    MeanBiexpLOOCVMinusKurtosisLOOCVArray(sf,d) = mean(CurBiexpLOOCVMinusKurtosisLOOCV(:));
    MeanBiexpLOOCVMinusMonoexpLOOCVArray(sf,d) = mean(CurBiexpLOOCVMinusMonoexpLOOCV(:));
    MeanKurtosisLOOCVMinusMonoexpLOOCVArray(sf,d) = mean(CurKurtosisLOOCVMinusMonoexpLOOCV(:));

    % Find the lowest value
    AICScores = [MonoexpAICArray(i,:); KurtAICArray(i,:); BiexpAICArray(i,:);];
    [~,LowestAICScores] = min(AICScores);
    % Monoexp = 1, Kurtosis = 2, Biexp = 3;
    LowestAICScorePctArray(sf,d,1) = sum(LowestAICScores == 1)/NoisySigDim * 100.; % Monoexp
    LowestAICScorePctArray(sf,d,2) = sum(LowestAICScores == 2)/NoisySigDim * 100.; % Kurtosis
    LowestAICScorePctArray(sf,d,3) = sum(LowestAICScores == 3)/NoisySigDim * 100.; % Biexp
    % AICc
    AICcScores = [MonoexpAICcArray(i,:); KurtAICcArray(i,:); BiexpAICcArray(i,:);];
    [~,LowestAICcScores] = min(AICcScores);
    % Monoexp = 1, Kurtosis = 2, Biexp = 3;
    LowestAICcScorePctArray(sf,d,1) = sum(LowestAICcScores == 1)/NoisySigDim * 100.; % Monoexp
    LowestAICcScorePctArray(sf,d,2) = sum(LowestAICcScores == 2)/NoisySigDim * 100.; % Kurtosis
    LowestAICcScorePctArray(sf,d,3) = sum(LowestAICcScores == 3)/NoisySigDim * 100.; % Biexp
    % LOOCV
    LOOCVScores = [MonoexpLOOCVScoreArray(i,:); KurtLOOCVScoreArray(i,:); BiexpLOOCVScoreArray(i,:);];
    [~,LowestLOOCVScores] = min(LOOCVScores);
    % Monoexp = 1, Kurtosis = 2, Biexp = 3;
    LowestLOOCVScorePctArray(sf,d,1) = sum(LowestLOOCVScores == 1)/NoisySigDim * 100.; % Monoexp
    LowestLOOCVScorePctArray(sf,d,2) = sum(LowestLOOCVScores == 2)/NoisySigDim * 100.; % Kurtosis
    LowestLOOCVScorePctArray(sf,d,3) = sum(LowestLOOCVScores == 3)/NoisySigDim * 100.; % Biexp
    
    % Head to head, biexp vs mono - AIC
    HeadToHeadAICBiexpVsMonoexp = [MonoexpAICArray(i,:); BiexpAICArray(i,:);];
    [~,LowestBiexpVsMonoexpAIC] = min(HeadToHeadAICBiexpVsMonoexp);
    BiexpAICSelectedVersusMonoexpArray(sf,d) = sum(LowestBiexpVsMonoexpAIC == 2)/NoisySigDim * 100;
    % and AICc
    HeadToHeadAICcBiexpVsMonoexp = [MonoexpAICcArray(i,:); BiexpAICcArray(i,:);];
    [~,LowestBiexpVsMonoexpAICc] = min(HeadToHeadAICcBiexpVsMonoexp);
    BiexpAICcSelectedVersusMonoexpArray(sf,d) = sum(LowestBiexpVsMonoexpAICc == 2)/NoisySigDim * 100;   
    

    % F-test results
    FNumerator = (MonoexpRSS11BVArray(i,:)-BiexpRSS11BVArray(i,:))./BiexpRSS11BVArray(i,:);
    FDenominator = (BiexpParam-MonoexpParam)./(BDim-BiexpParam);
    FRatio = FNumerator./FDenominator;

    FCriticalValue = finv(0.95, BiexpParam-MonoexpParam, BDim-BiexpParam);
    FtestBiexpBetter = find(FRatio > FCriticalValue);
%     FTestScoreArray = zeros(SF1Dim, D1D2Dim);
    FTestBiexpBetterPctArray(sf,d) = length(FtestBiexpBetter)/NoisySigDim * 100;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hard coded from the parameter creation code
SF1Values = 0:1/69:1;
D1D2RatioValues = 2:18/69:20;
PlotFontSize = 12;
% Compare AIC and F-Test
figure('Position', [0,0,1000,800]);
subplot(2,2,1);h=pcolor(SF1Values, D1D2RatioValues, BiexpAICSelectedVersusMonoexpArray(:,:)');set(h,'edgecolor','none');title('Biexp AIC Selection % vs. Monoexp','fontsize', PlotFontSize+1);hold on;
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
HandleC = colorbar('location', 'EastOutside'); HandleCLabel = ylabel(HandleC,'%'); drawnow; % MUST HAVE
LabelPos = HandleCLabel.Position; HandleCLabel.Rotation = -90; HandleCLabel.Position = [LabelPos(1) + 0.4, LabelPos(2), LabelPos(3)];

subplot(2,2,2);h=pcolor(SF1Values, D1D2RatioValues, FTestBiexpBetterPctArray(:,:)');set(h,'edgecolor','none');title('Biexp \itF\rm\bf-Test Selection % vs Monoexp','fontsize', PlotFontSize+1);hold on;
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
HandleC = colorbar('location', 'EastOutside'); HandleCLabel = ylabel(HandleC,'%'); drawnow; % MUST HAVE
LabelPos = HandleCLabel.Position; HandleCLabel.Rotation = -90; HandleCLabel.Position = [LabelPos(1) + 0.4, LabelPos(2), LabelPos(3)];

subplot(2,2,4);h=pcolor(SF1Values, D1D2RatioValues, BiexpAICcSelectedVersusMonoexpArray(:,:)');set(h,'edgecolor','none');title('Biexp AIC_{c} Selection % vs. Monoexp','fontsize', PlotFontSize+1);hold on;set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
HandleC = colorbar('location', 'EastOutside'); HandleCLabel = ylabel(HandleC,'%'); drawnow; % MUST HAVE
LabelPos = HandleCLabel.Position; HandleCLabel.Rotation = -90; HandleCLabel.Position = [LabelPos(1) + 0.4, LabelPos(2), LabelPos(3)];
