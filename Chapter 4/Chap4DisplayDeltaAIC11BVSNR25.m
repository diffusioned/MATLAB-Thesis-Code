% Chap4DisplayDeltaAIC11BVSNR25.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap4DisplayDeltaAIC11BVSNR25()

% Load noise free parameter values for data organization
load('YourPath\BiexpNoiseFreeParameters.mat');

% Load model fits
load('YourPath\BiexpFitArray11BValues.mat');
load('YourPath\MonoexpFitArray11BValues.mat');
load('YourPath\KurtFitArray11BValues.mat');
PlotFontSize = 10;

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
LowestAICScorePctArray  = zeros(SF1Dim, D1D2Dim, 3, 2); % For 3 model combinations

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

    % Separate Test Selection Percentages by two model combinations
    AICScoresBiexpKurt = [KurtAICArray(i,:); BiexpAICArray(i,:);];
    AICScoresBiexpMonoexp = [MonoexpAICArray(i,:); BiexpAICArray(i,:);];
    AICScoresKurtMonoexp = [MonoexpAICArray(i,:); KurtAICArray(i,:);];
    [~,LowestAICScoresBiexpKurt] = min(AICScoresBiexpKurt);
    [~,LowestAICScoresBiexpMonoexp] = min(AICScoresBiexpMonoexp);
    [~,LowestAICScoresKurtMonoexp] = min(AICScoresKurtMonoexp);
    % Monoexp = 1, Kurtosis = 2, Biexp = 3;
    LowestAICScorePctArray(sf,d,1,1) = sum(LowestAICScoresBiexpKurt == 1)/NoisySigDim * 100.; % Kurt
    LowestAICScorePctArray(sf,d,1,2) = sum(LowestAICScoresBiexpKurt == 2)/NoisySigDim * 100.; % Biexp
    
    LowestAICScorePctArray(sf,d,2,1) = sum(LowestAICScoresBiexpMonoexp == 1)/NoisySigDim * 100.; % Monoexp
    LowestAICScorePctArray(sf,d,2,2) = sum(LowestAICScoresBiexpMonoexp == 2)/NoisySigDim * 100.; % Biexp
    
    LowestAICScorePctArray(sf,d,3,1) = sum(LowestAICScoresKurtMonoexp == 1)/NoisySigDim * 100.; % Monoexp
    LowestAICScorePctArray(sf,d,3,2) = sum(LowestAICScoresKurtMonoexp == 2)/NoisySigDim * 100.; % Kurt
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SF1Values = 0:1/69:1;
D1D2RatioValues = 2:18/69:20;
SideBySidePlotFontSize = 8;
% AIC Biexp - Kurtosis
figure('Position', [50,50,700,800]);subplot(3,2,1);
h=pcolor(SF1Values, D1D2RatioValues, MeanBiexpAICMinusKurtosisAICArray(:,:)');colorbar;set(h,'edgecolor','none');title('\it\DeltaAIC(Biexp - Kurtosis)', 'fontsize', SideBySidePlotFontSize+2);hold on;
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
contour(SF1Values, D1D2RatioValues, MeanBiexpAICMinusKurtosisAICArray(:,:)', [0,0], 'LineColor', [0 0 0], 'LineWidth', 1.5);
% AIC Biexp - Monoexp
subplot(3,2,3);h=pcolor(SF1Values, D1D2RatioValues, MeanBiexpAICMinusMonoexpAICArray(:,:)');colorbar;set(h,'edgecolor','none');title('\it\DeltaAIC(Biexp - Mono)', 'fontsize', SideBySidePlotFontSize+2);hold on;
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
contour(SF1Values, D1D2RatioValues, MeanBiexpAICMinusMonoexpAICArray(:,:)', [0,0], 'LineColor', [0 0 0], 'LineWidth', 1.5);
% AIC Kurtosis - Monoexp
subplot(3,2,5);h=pcolor(SF1Values, D1D2RatioValues, MeanKurtosisAICMinusMonoexpAICArray(:,:)');colorbar;set(h,'edgecolor','none');title('\it\DeltaAIC(Kurtosis - Mono)', 'fontsize', SideBySidePlotFontSize+2);hold on;
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
contour(SF1Values, D1D2RatioValues, MeanKurtosisAICMinusMonoexpAICArray(:,:)', [0,0], 'LineColor', [0 0 0], 'LineWidth', 1.5);

% Selection percentages
% AIC Selection Pct for Biexp vs Kurt
subplot(3,2,2);h=pcolor(SF1Values, D1D2RatioValues, squeeze(LowestAICScorePctArray(:,:,1,2))');colorbar;set(h,'edgecolor','none');title('Biexp AIC Selection % vs Kurt', 'fontsize', SideBySidePlotFontSize+2);caxis([0 100]);hold on;
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
contour(SF1Values, D1D2RatioValues, LowestAICScorePctArray(:,:,1,2)', [50,50], 'LineColor', [0 0 0], 'LineWidth', 1.5);
% AIC Selection Pct for Biexp vs Mono
subplot(3,2,4);h=pcolor(SF1Values, D1D2RatioValues, squeeze(LowestAICScorePctArray(:,:,2,2))');colorbar;set(h,'edgecolor','none');title('Biexp AIC Selection % vs Mono', 'fontsize', SideBySidePlotFontSize+2);caxis([0 100]);hold on;
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
contour(SF1Values, D1D2RatioValues, LowestAICScorePctArray(:,:,2,2)', [50,50], 'LineColor', [0 0 0], 'LineWidth', 1.5);
% AIC Selection Pct for Kurt vs Mono
subplot(3,2,6);h=pcolor(SF1Values, D1D2RatioValues, squeeze(LowestAICScorePctArray(:,:,3,2))');colorbar;set(h,'edgecolor','none');title('Kurt AIC Selection % vs Mono', 'fontsize', SideBySidePlotFontSize+2);caxis([0 100]);hold on;
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
contour(SF1Values, D1D2RatioValues, LowestAICScorePctArray(:,:,3,2)', [50,50], 'LineColor', [0 0 0], 'LineWidth', 1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare delAIC and delAICc
% AIC Biexp - Kurtosis
SideBySidePlotFontSize = 8;
figure('Position', [0,0,700,800]);subplot(3,2,1);
h=pcolor(SF1Values, D1D2RatioValues, MeanBiexpAICMinusKurtosisAICArray(:,:)');colorbar;set(h,'edgecolor','none');title('\it\DeltaAIC(Biexp - Kurtosis)', 'fontsize', SideBySidePlotFontSize+2);hold on;
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
contour(SF1Values, D1D2RatioValues, MeanBiexpAICMinusKurtosisAICArray(:,:)', [0,0], 'LineColor', [0 0 0], 'LineWidth', 1.5);
% AIC Biexp - Monoexp
subplot(3,2,3);h=pcolor(SF1Values, D1D2RatioValues, MeanBiexpAICMinusMonoexpAICArray(:,:)');colorbar;set(h,'edgecolor','none');title('\it\DeltaAIC(Biexp -  Monoexp)', 'fontsize', SideBySidePlotFontSize+2);hold on;
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
contour(SF1Values, D1D2RatioValues, MeanBiexpAICMinusMonoexpAICArray(:,:)', [0,0], 'LineColor', [0 0 0], 'LineWidth', 1.5);
% AIC Kurtosis - Monoexp
subplot(3,2,5);h=pcolor(SF1Values, D1D2RatioValues, MeanKurtosisAICMinusMonoexpAICArray(:,:)');colorbar;set(h,'edgecolor','none');title('\it\DeltaAIC(Kurtosis -  Monoexp)', 'fontsize', SideBySidePlotFontSize+2);hold on;
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
contour(SF1Values, D1D2RatioValues, MeanBiexpAICMinusMonoexpAICArray(:,:)', [0,0], 'LineColor', [0 0 0], 'LineWidth', 1.5);

% AICc Biexp - Kurtosis
subplot(3,2,2);h=pcolor(SF1Values, D1D2RatioValues, MeanBiexpAICcMinusKurtosisAICcArray(:,:)');colorbar;set(h,'edgecolor','none');title('\it\DeltaAIC_{c}(Biexp - Kurtosis)', 'fontsize', SideBySidePlotFontSize+2);hold on;
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
contour(SF1Values, D1D2RatioValues, MeanBiexpAICcMinusKurtosisAICcArray(:,:)', [0,0], 'LineColor', [0 0 0], 'LineWidth', 1.5);
% AICc Biexp - Monoexp
subplot(3,2,4);h=pcolor(SF1Values, D1D2RatioValues, MeanBiexpAICcMinusMonoexpAICcArray(:,:)');colorbar;set(h,'edgecolor','none');title('\it\DeltaAIC_{c}(Biexp -  Monoexp)', 'fontsize', SideBySidePlotFontSize+2);hold on;
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
contour(SF1Values, D1D2RatioValues, MeanBiexpAICcMinusMonoexpAICcArray(:,:)', [0,0], 'LineColor', [0 0 0], 'LineWidth', 1.5);
% AICc Kurtosis - Monoexp
subplot(3,2,6);h=pcolor(SF1Values, D1D2RatioValues, MeanKurtosisAICcMinusMonoexpAICcArray(:,:)');colorbar;set(h,'edgecolor','none');title('\it\DeltaAIC_{c}(Kurtosis -  Monoexp)', 'fontsize', SideBySidePlotFontSize+2);hold on;
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', SideBySidePlotFontSize);
contour(SF1Values, D1D2RatioValues, MeanKurtosisAICcMinusMonoexpAICcArray(:,:)', [0,0], 'LineColor', [0 0 0], 'LineWidth', 1.5);
