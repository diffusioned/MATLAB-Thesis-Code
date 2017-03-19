% Chap4DisplayDeltaAICCaseStudies11BVSNR25.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap4DisplayDeltaAICCaseStudies11BVSNR25()

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
    
    % Head to head
    HeadToHeadAICBiexpVsMonoexp = [MonoexpAICArray(i,:); BiexpAICArray(i,:);];
    [~,LowestBiexpVsMonoexpAIC] = min(HeadToHeadAICBiexpVsMonoexp);
    BiexpAICSelectedVersusMonoexpArray(sf,d) = sum(LowestBiexpVsMonoexpAIC == 2)/NoisySigDim * 100;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the biexp vs monoexp selection rate and compare the difference in delAIC values for 3 cases
% Hard coded from the parameter creation code
SF1Values = 0:1/69:1;
D1D2RatioValues = 2:18/69:20;
% Biexp vs Monoexp
figure('Position', [0,0,1000,800]);
subplot(2,2,1);h=pcolor(SF1Values, D1D2RatioValues, BiexpAICSelectedVersusMonoexpArray(:,:)');colorbar;set(h,'edgecolor','none');title('Biexp AIC Selection % vs. Mono');hold on;
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% contour(SF1Values, D1D2RatioValues, MeanBiexpAICMinusMonoexpAICArray(:,:)', [0,0], 'LineColor', [1 1 1], 'LineWidth', 1);
text(0.4928,19.7,'A','Color','r','FontSize',18, 'FontWeight', 'bold'); % Signal Index 4865
text(0.087,19.1,'B','Color','r','FontSize',18, 'FontWeight', 'bold'); % Signal Index 4697
text(-0.02,13.22,'C','Color','r','FontSize',18, 'FontWeight', 'bold'); % Signal Index 4865

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note - your values may differ for your set, so change the specific numbers below for the three signals
% SIGNAL A
% 4865 = 100% selection of biexponential, distinct difference between distributions, SF1Idx = 35, D1D2Idx = 70
CurBiexpMinusMonoexpAIC4865 = squeeze(BiexpAICMinusMonoexpAICArray(4865,:));
% The distribution of the difference here is near normal
subplot(2,2,2);histogram(BiexpAICMinusMonoexpAICArray(4865,:),30);
title('{\color{red}A}: 100% Biexp Selection Rate','FontWeight','bold');
xlabel('\it\DeltaAIC(Biexp - Mono)','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Count','FontWeight','bold','fontsize', PlotFontSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGNAL B
% Specific example with selection rate close to 50%  SF1Idx = 7, D1D2Idx = 68
NFSignalIdx = 4697; % selection rate = 50%, mean diff AIC = -0.9449
CurBiexpMinusMonoexpAI4697 = squeeze(BiexpAICMinusMonoexpAICArray(4697,:));
subplot(2,2,3);histogram(BiexpAICMinusMonoexpAICArray(4697,:),30);
title('{\color{red}B}: 50% Biexp Selection Rate');
xlabel('\it\DeltaAIC(Biexp - Mono)','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Count','FontWeight','bold','fontsize', PlotFontSize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGNAL C
% Indices:
% 3011 = 5% selection of biexponential, distinct difference between distributions, SF1Idx = 1, D1D2Idx = 44
CurBiexpMinusMonoexpAIC3011 = squeeze(BiexpAICMinusMonoexpAICArray(3011,:));
subplot(2,2,4);histogram(BiexpAICMinusMonoexpAICArray(3011,:),30);
title('{\color{red}C}: 5% Biexp Selection Rate');
xlabel('\it\DeltaAIC(Biexp - Mono)','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Count','FontWeight','bold','fontsize', PlotFontSize);

% Mean values for reference
mean(BiexpAICMinusMonoexpAICArray(4865,:))
mean(BiexpAICMinusMonoexpAICArray(4697,:))
mean(BiexpAICMinusMonoexpAICArray(3011,:))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do the same 3 cases, only displaying a 3x2 image of the RSS and AIC values of the biexponential and monoexponental fits
figure('Position', [0,0,800,800]);
% SIGNAL A
% 4865 = 100% selection of biexponential, distinct difference between distributions, SF1Idx = 35, D1D2Idx = 70
subplot(3,2,1);histogram(BiexpAICArray(4865,:),-90:1:-30);
hold all;histogram(MonoexpAICArray(4865,:),-90:1:-30);
title('Signal A AIC: {\color{blue}Biexp} vs. {\color{red}Mono}','FontWeight','bold');
xlabel('AIC','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Count','FontWeight','bold','fontsize', PlotFontSize);

subplot(3,2,2);histogram(BiexpRSS11BVArray(4865,:),0:0.005:0.25);
hold all;histogram(MonoexpRSS11BVArray(4865,:),0:0.005:0.25);
title('Signal A RSS: {\color{blue}Biexp} vs. {\color{red}Mono}','FontWeight','bold');
xlabel('RSS','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Count','FontWeight','bold','fontsize', PlotFontSize);

subplot(3,2,3);histogram(BiexpAICArray(4697,:),-90:0.75:-50);
hold all;histogram(MonoexpAICArray(4697,:),-90:0.75:-50);
title('Signal B AIC: {\color{blue}Biexp} vs. {\color{red}Mono}','FontWeight','bold');
xlabel('AIC','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Count','FontWeight','bold','fontsize', PlotFontSize);

subplot(3,2,4);histogram(BiexpRSS11BVArray(4697,:),0:0.001:0.05);
hold all;histogram(MonoexpRSS11BVArray(4697,:),0:0.001:0.025);
title('Signal B RSS: {\color{blue}Biexp} vs. {\color{red}Mono}','FontWeight','bold');
xlabel('RSS','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Count','FontWeight','bold','fontsize', PlotFontSize);

subplot(3,2,5);histogram(BiexpAICArray(3011,:),-100:1:-50);
hold all;histogram(MonoexpAICArray(3011,:),-100:1:-50);
title('Signal C AIC: {\color{blue}Biexp} vs. {\color{red}Mono}','FontWeight','bold');
xlabel('AIC','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Count','FontWeight','bold','fontsize', PlotFontSize);

subplot(3,2,6);histogram(BiexpRSS11BVArray(3011,:),0:0.001:0.05);
hold all;histogram(MonoexpRSS11BVArray(3011,:),0:0.001:0.05);
title('Signal C RSS: {\color{blue}Biexp} vs. {\color{red}Mono}','FontWeight','bold');
xlabel('RSS','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Count','FontWeight','bold','fontsize', PlotFontSize);