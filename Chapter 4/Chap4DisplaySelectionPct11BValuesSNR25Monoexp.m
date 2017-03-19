% Chap4DisplaySelectionPct11BValuesSNR25Monoexp.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap4DisplaySelectionPct11BValuesSNR25Monoexp()

% Load noise free parameter values for data organization
load('YourPath\MonoexpNFParametersForKurtosis.mat');

% Load model fits
load('YourPath\BiexpFitArrayMONO11BValues.mat');
load('YourPath\MonoexpFitArrayMONO11BValues.mat');
load('YourPath\KurtFitArrayMONO11BValues.mat');
PlotFontSize = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load diagnostic metrics for monoexp fits
load('YourPath\BiexpRSSMONO11BVArray.mat');
load('YourPath\BiexpLOOCVScoreArrayMONO11BValues.mat');

% Kurtosis fits
load('YourPath\KurtRSSMONO11BVArray.mat');
load('YourPath\KurtLOOCVScoreArrayMONO11BValues.mat');

% Monoexp fits
load('YourPath\MonoexpRSSMONO11BVArray.mat');
load('YourPath\MonoexpLOOCVScoreArrayMONO11BValues.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare the AIC scores for monoexp truth and display the plots
% Grab the dimensions
[NoiseFreeDim, ~] = size(NFTestParameterArray);
[~, NoisySigDim, ~] = size(BiexpFitArray);

% Explicitly set
BDim = 11;

% Sort the noise free parameters into bins - create 96 bins
ADCBinEdges = 0.05:0.01:1.; ADCBinDim = 95;
EdgeVals = {ADCBinEdges}; % Cell array
[ADCCnt,~,ADCIdxs] = histcounts(NFTestParameterArray(:,2),ADCBinEdges);

% Calculate the AIC values
MonoexpParam = 2+1; KurtParam = 3+1;  BiexpParam = 4+1; % The +1 is to include the variance needed with RSS fits
MonoexpAICArray = BDim.*log(MonoexpRSS11BVArray/BDim) + 2*(MonoexpParam);
KurtAICArray = BDim.*log(KurtRSS11BVArray/BDim) + 2*(KurtParam);
BiexpAICArray = BDim.*log(BiexpRSS11BVArray/BDim) + 2*(BiexpParam);

% Create arrays for statistical calculations
LowestAICScorePctArray  = zeros(ADCBinDim, 3); % For 3 models

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AICc - repeat for corrected AIC values
% Calculate the AICc values
MonoexpAICcArray = MonoexpAICArray + 2*MonoexpParam*(MonoexpParam+1)/(BDim-MonoexpParam-1); % 11 B values = 3.4286
KurtAICcArray = KurtAICArray + 2*KurtParam*(KurtParam+1)/(BDim-KurtParam-1); % 11 B values = 6.6667
BiexpAICcArray = BiexpAICArray + 2*BiexpParam*(BiexpParam+1)/(BDim-BiexpParam-1); % 11 B values = 12
LowestAICcScorePctArray  = zeros(ADCBinDim, 3); % For 3 models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assessment of LOOCV scores
LowestLOOCVScorePctArray  = zeros(ADCBinDim, 3); % For 3 models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop and calculate statisticals for each noise-free signal
for i = 1:ADCBinDim
    % Select the noise frees signals in this bin
    CurBinIdxs = find(ADCIdxs == i);
    CurMonoexpAIC = MonoexpAICArray(CurBinIdxs,:);  CurKurtAIC = KurtAICArray(CurBinIdxs,:);  CurBiexpAIC = BiexpAICArray(CurBinIdxs,:);
    CurMonoexpAICc = MonoexpAICcArray(CurBinIdxs,:);  CurKurtAICc = KurtAICcArray(CurBinIdxs,:);  CurBiexpAICc = BiexpAICcArray(CurBinIdxs,:);
    CurMonoexpLOOCV = MonoexpLOOCVScoreArray(CurBinIdxs,:);  CurKurtLOOCV = KurtLOOCVScoreArray(CurBinIdxs,:);  CurBiexpLOOCV = BiexpLOOCVScoreArray(CurBinIdxs,:);
    
    % Find the lowest value.  Note ' below to change the array shape
    AICScores = [CurMonoexpAIC(:)'; CurKurtAIC(:)'; CurBiexpAIC(:)';];
    [~,LowestAICScores] = min(AICScores);
    BinSigDim = length(CurBinIdxs) * NoisySigDim;
    % Monoexp = 1, Kurtosis = 2, Biexp = 3;
    LowestAICScorePctArray(i,1) = sum(LowestAICScores == 1)/BinSigDim * 100; % Monoexp
    LowestAICScorePctArray(i,2) = sum(LowestAICScores == 2)/BinSigDim * 100; % Kurtosis
    LowestAICScorePctArray(i,3) = sum(LowestAICScores == 3)/BinSigDim * 100; % Biexp
    % AICc
    AICcScores = [CurMonoexpAICc(:)'; CurKurtAICc(:)'; CurBiexpAICc(:)';];
    [~,LowestAICcScores] = min(AICcScores);
    % Monoexp = 1, Kurtosis = 2, Biexp = 3;
    LowestAICcScorePctArray(i,1) = sum(LowestAICcScores == 1)/BinSigDim * 100; % Monoexp
    LowestAICcScorePctArray(i,2) = sum(LowestAICcScores == 2)/BinSigDim * 100; % Kurtosis
    LowestAICcScorePctArray(i,3) = sum(LowestAICcScores == 3)/BinSigDim * 100; % Biexp
    % LOOCV
    LOOCVScores = [CurMonoexpLOOCV(:)'; CurKurtLOOCV(:)'; CurBiexpLOOCV(:)';];
    [~,LowestLOOCVScores] = min(LOOCVScores);
    % Monoexp = 1, Kurtosis = 2, Biexp = 3;
    LowestLOOCVScorePctArray(i,1) = sum(LowestLOOCVScores == 1)/BinSigDim * 100; % Monoexp
    LowestLOOCVScorePctArray(i,2) = sum(LowestLOOCVScores == 2)/BinSigDim * 100; % Kurtosis
    LowestLOOCVScorePctArray(i,3) = sum(LowestLOOCVScores == 3)/BinSigDim * 100; % Biexp
    i
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single maps for thesis images 11 B Values
SF1Values = 0:1/69:1; % LEAVE FOR ALL BELOW
D1D2RatioValues = 2:18/69:20; % LEAVE FOR ALL BELOW
PlotFontSize = 18; % LEAVE FOR ALL BELOW

% AIC Selection Pct for Monoexp
figure('Position', [0,0,800,800]);
h1=bar(ADCBinEdges(2:end),squeeze(LowestAICScorePctArray(:,1)),'histc');
xlim([0.06 1.0]);ylim([0 100]);
xlabel('True \itADC\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Selection Percentage','FontWeight','bold','fontsize', PlotFontSize);
% title('Monoexponential Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

% AIC Selection Pct for Kurtosis
figure('Position', [0,0,800,800]);
h2=bar(ADCBinEdges(2:end),squeeze(LowestAICScorePctArray(:,2)),'histc');
xlim([0.06 1.0]);ylim([0 100]);
xlabel('True \itADC\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Selection Percentage','FontWeight','bold','fontsize', PlotFontSize);
% title('Kurtosis Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

% AIC Selection Pct for Biexp
figure('Position', [0,0,800,800]);
h3=bar(ADCBinEdges(2:end),squeeze(LowestAICScorePctArray(:,3)),'histc');
xlim([0.06 1.0]);ylim([0 100]);
xlabel('True \itADC\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Selection Percentage','FontWeight','bold','fontsize', PlotFontSize);
% title('Biexponential Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

%%%%% AICc %%%%%
% AICc Selection Pct for Monoexp
figure('Position', [0,0,800,800]);
h4=bar(ADCBinEdges(2:end),squeeze(LowestAICcScorePctArray(:,1)),'histc');
xlim([0.06 1.0]);ylim([0 100]);
xlabel('True \itADC\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Selection Percentage','FontWeight','bold','fontsize', PlotFontSize);
% title('Monoexponential Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

% AICc Selection Pct for Kurtosis
figure('Position', [0,0,800,800]);
h5=bar(ADCBinEdges(2:end),squeeze(LowestAICcScorePctArray(:,2)),'histc');
xlim([0.06 1.0]);ylim([0 100]);
xlabel('True \itADC\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Selection Percentage','FontWeight','bold','fontsize', PlotFontSize);
% title('Kurtosis Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

% AICc Selection Pct for Biexp
figure('Position', [0,0,800,800]);
h6=bar(ADCBinEdges(2:end),squeeze(LowestAICcScorePctArray(:,3)),'histc');
xlim([0.06 1.0]);ylim([0 100]);
xlabel('True \itADC\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Selection Percentage','FontWeight','bold','fontsize', PlotFontSize);
% title('Biexponential Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

%%%%% LOOCV %%%%%
% LOOCV Selection Pct for Monoexp
figure('Position', [0,0,800,800]);
h7=bar(ADCBinEdges(2:end),squeeze(LowestLOOCVScorePctArray(:,1)),'histc');
xlim([0.06 1.0]);ylim([0 100]);
xlabel('True \itADC\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Selection Percentage','FontWeight','bold','fontsize', PlotFontSize);
% title('Monoexponential Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

% LOOCV Selection Pct for Kurtosis
figure('Position', [0,0,800,800]);
h8=bar(ADCBinEdges(2:end),squeeze(LowestLOOCVScorePctArray(:,2)),'histc');
xlim([0.06 1.0]);ylim([0 100]);
xlabel('True \itADC\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Selection Percentage','FontWeight','bold','fontsize', PlotFontSize);
% title('Kurtosis Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

% LOOCV Selection Pct for Biexp
figure('Position', [0,0,800,800]);
h9=bar(ADCBinEdges(2:end),squeeze(LowestLOOCVScorePctArray(:,3)),'histc');
xlim([0.06 1.0]);ylim([0 100]);
xlabel('True \itADC\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('Selection Percentage','FontWeight','bold','fontsize', PlotFontSize);
% title('Biexponential Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
