% Chap4DisplaySelectionPct11BValuesSNR25Biexp.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap4DisplaySelectionPct11BValuesSNR25Biexp()

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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single maps for thesis images 11 B Values
SF1Values = 0:1/69:1; % LEAVE FOR ALL BELOW
D1D2RatioValues = 2:18/69:20; % LEAVE FOR ALL BELOW
PlotFontSize = 18; % LEAVE FOR ALL BELOW

% AIC Selection Pct for Monoexp
figure('Position', [0,0,800,800]);
h=pcolor(SF1Values, D1D2RatioValues, squeeze(LowestAICScorePctArray(:,:,1))');set(h,'edgecolor','none');
caxis([0 100]);
HandleC =  colorbar('location', 'NorthOutside');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% title('Monoexponential Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

% AIC Selection Pct for Kurtosis
figure('Position', [0,0,800,800]);
h=pcolor(SF1Values, D1D2RatioValues, squeeze(LowestAICScorePctArray(:,:,2))');set(h,'edgecolor','none');
caxis([0 100]);
HandleC =  colorbar('location', 'NorthOutside');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% title('Kurtosis Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

% AIC Selection Pct for Biexp
figure('Position', [0,0,800,800]);
h=pcolor(SF1Values, D1D2RatioValues, squeeze(LowestAICScorePctArray(:,:,3))');set(h,'edgecolor','none');
caxis([0 100]);
HandleC =  colorbar('location', 'NorthOutside');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% title('Biexponential Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

%%%%% AICc %%%%%
% AICc Selection Pct for Monoexp
figure('Position', [0,0,800,800]);
h=pcolor(SF1Values, D1D2RatioValues, squeeze(LowestAICcScorePctArray(:,:,1))');set(h,'edgecolor','none');
caxis([0 100]);
HandleC =  colorbar('location', 'NorthOutside');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% title('Monoexponential Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

% AICc Selection Pct for Kurtosis
figure('Position', [0,0,800,800]);
h=pcolor(SF1Values, D1D2RatioValues, squeeze(LowestAICcScorePctArray(:,:,2))');set(h,'edgecolor','none');
caxis([0 100]);
HandleC =  colorbar('location', 'NorthOutside');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% title('Kurtosis Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

% AICc Selection Pct for Biexp
figure('Position', [0,0,800,800]);
h=pcolor(SF1Values, D1D2RatioValues, squeeze(LowestAICcScorePctArray(:,:,3))');set(h,'edgecolor','none');
caxis([0 100]);
HandleC =  colorbar('location', 'NorthOutside');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% title('Biexponential Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

%%%%% LOOCV %%%%%
% LOOCV Selection Pct for Monoexp
figure('Position', [0,0,800,800]);
h=pcolor(SF1Values, D1D2RatioValues, squeeze(LowestLOOCVScorePctArray(:,:,1))');set(h,'edgecolor','none');
caxis([0 100]);
HandleC =  colorbar('location', 'NorthOutside');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% title('Monoexponential Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

% LOOCV Selection Pct for Kurtosis
figure('Position', [0,0,800,800]);
h=pcolor(SF1Values, D1D2RatioValues, squeeze(LowestLOOCVScorePctArray(:,:,2))');set(h,'edgecolor','none');
caxis([0 100]);
HandleC =  colorbar('location', 'NorthOutside');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itADC_{1} \rm\bf/ \itADC_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% title('Kurtosis Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size

% LOOCV Selection Pct for Biexp
figure('Position', [0,0,800,800]);
h=pcolor(SF1Values, D1D2RatioValues, squeeze(LowestLOOCVScorePctArray(:,:,3))');set(h,'edgecolor','none');
caxis([0 100]);
HandleC =  colorbar('location', 'NorthOutside');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itADC_{1} \rm\bf/ \itADC_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% title('Biexponential Selection Percentage');
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
