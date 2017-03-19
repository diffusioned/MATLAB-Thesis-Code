% Chap2PlotHistograms.m
% MATLAB file for displaying histograms in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap2PlotHistograms()
% A little explanation here.  I first loaded my noise free parameter set to pick out values
load('YourPath\YourNoiseFreeParameters.mat');
% First test set - D1/D2 = 15, SF1 = 0.05, 0.1, 0.25, 0.5

% Use the find code to fill in the index numbers
% Closest to SF = 0.5
find(NFTestParameterArray(:,1) > 0.501 & NFTestParameterArray(:,1) < 0.502 & 1./NFTestParameterArray(:,4) > 14.9 & 1./NFTestParameterArray(:,4) < 15.1 )
SF05_D1D215_Index = 0000;
% Closest to SF = 0.25
find(NFTestParameterArray(:,1) > 0.249 & NFTestParameterArray(:,1) < 0.251 & 1./NFTestParameterArray(:,4) > 14.9 & 1./NFTestParameterArray(:,4) < 15.1 )
SF025_D1D215_Index = 0000;
% Closest to SF = 0.1
find(NFTestParameterArray(:,1) > 0.095 & NFTestParameterArray(:,1) < 0.105 & 1./NFTestParameterArray(:,4) > 14.95 & 1./NFTestParameterArray(:,4) < 15.05 )
SF01_D1D215_Index = 0000;
% Closest to SF = 0.05
find(NFTestParameterArray(:,1) > 0.045 & NFTestParameterArray(:,1) < 0.055 & 1./NFTestParameterArray(:,4) > 14.95 & 1./NFTestParameterArray(:,4) < 15.05 )
SF005_D1D215_Index = 0000;

% Load the corresponding parameter fits
load('YourPath\BiexpFitArray.mat');
% Load the covariance matrix to extract the parameter standard errors, if needed
load('YourPath\CovMxArray.mat');

% 3x4 Histograms plus added true value lines
figure('Position', [50,00,1000,1000]);
h1 = subplot(4,3,1);histogram(BiexpFitArray(SF05_D1D215_Index,:,1),0:0.02:1);title('\itSF_{1}\rm\bf = 0.5, \itD_{1} / D_{2}\rm\bf = 15');ylabel('\itA_{1}\rm\bf Estimates','FontWeight','bold');ylim([0 40]);xlim([0 1]);
hold on; plot(repmat(0.5,1,2),[0 40], 'linewidth',2, 'linestyle','--'); % or 1.5*max(h.Values) if needed for adjustable
h2 = subplot(4,3,4);histogram(BiexpFitArray(SF05_D1D215_Index,:,2),0:0.02:1);ylabel('\itA_{2}\rm\bf Estimates','FontWeight','bold');ylim([0 40]);xlim([0 1.5]);
hold on; plot(repmat(0.5,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h3 = subplot(4,3,7);histogram(BiexpFitArray(SF05_D1D215_Index,:,3),0:0.08:4);ylabel('\itD_{1}\rm\bf Estimates','FontWeight','bold');ylim([0 40]);xlim([0 4]);
hold on; plot(repmat(1,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h4 = subplot(4,3,10);histogram(BiexpFitArray(SF05_D1D215_Index,:,4),0:0.004:0.2);ylabel('\itD_{2}\rm\bf Estimates','FontWeight','bold');ylim([0 40]);xlim([-0.05 0.15]);
hold on; plot(repmat(0.067,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h5 = subplot(4,3,2);histogram(BiexpFitArray(SF025_D1D215_Index,:,1),0:0.02:1);title('\itSF_{1}\rm\bf = 0.25, \itD_{1} / D_{2}\rm\bf = 15');ylim([0 40]);xlim([0 1]);
hold on; plot(repmat(0.25,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h6 = subplot(4,3,5);histogram(BiexpFitArray(SF025_D1D215_Index,:,2),0:0.02:1);ylim([0 40]); xlim([0 1.5]);
hold on; plot(repmat(0.75,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h7 = subplot(4,3,8);histogram(BiexpFitArray(SF025_D1D215_Index,:,3),0:0.08:4);ylim([0 40]);xlim([0 4]);
hold on; plot(repmat(1,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h8 = subplot(4,3,11);histogram(BiexpFitArray(SF025_D1D215_Index,:,4),0:0.004:0.2);ylim([0 40]);xlim([-0.05 0.15]);
hold on; plot(repmat(0.067,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h9 = subplot(4,3,3);histogram(BiexpFitArray(SF01_D1D215_Index,:,1),0:0.02:1);title('\itSF_{1}\rm\bf = 0.1, \itD_{1} / D_{2}\rm\bf = 15');ylim([0 40]);xlim([-0.1 1]);
hold on; plot(repmat(0.1,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h10 = subplot(4,3,6);histogram(BiexpFitArray(SF01_D1D215_Index,:,2),0:0.02:1);ylim([0 40]);xlim([0 1.5]);
hold on; plot(repmat(0.9,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h11 = subplot(4,3,9);histogram(BiexpFitArray(SF01_D1D215_Index,:,3),0:0.08:4);ylim([0 40]);xlim([-0.1 4]);
hold on; plot(repmat(1,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h12 = subplot(4,3,12);histogram(BiexpFitArray(SF01_D1D215_Index,:,4),0:0.004:0.2);ylim([0 40]); xlim([-0.5 0.5]);
hold on; plot(repmat(0.067,1,2),[0 40], 'linewidth',2, 'linestyle','--');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD CONFIDENCE INTERVALS IF NEEDED
DOFModel = 11-4; % signals - fit parameters
alpha = 0.05; % 95% CI
TFactor = tinv(1-alpha/2,DOFModel);

CurIdxs = [SF05_D1D215_Index,SF025_D1D215_Index,SF01_D1D215_Index];
CurSEAmp1Array = sqrt(CovMxArray(CurIdxs,:,1,1));
CurSEAmp2Array = sqrt(CovMxArray(CurIdxs,:,2,2));
CurSEADC1Array = sqrt(CovMxArray(CurIdxs,:,3,3));
CurSEADC2Array = sqrt(CovMxArray(CurIdxs,:,4,4));

MedianMinAmp1CIValues = median(squeeze(BiexpFitArray(CurIdxs,:,1)) - CurSEAmp1Array*TFactor,2);
MedianMinAmp2CIValues = median(squeeze(BiexpFitArray(CurIdxs,:,2)) - CurSEAmp2Array*TFactor,2);
MedianMinADC1CIValues = median(squeeze(BiexpFitArray(CurIdxs,:,3)) - CurSEADC1Array*TFactor,2);
MedianMinADC2CIValues = median(squeeze(BiexpFitArray(CurIdxs,:,4)) - CurSEADC2Array*TFactor,2);
MedianMaxAmp1CIValues = median(squeeze(BiexpFitArray(CurIdxs,:,1)) + CurSEAmp1Array*TFactor,2);
MedianMaxAmp2CIValues = median(squeeze(BiexpFitArray(CurIdxs,:,2)) + CurSEAmp2Array*TFactor,2);
MedianMaxADC1CIValues = median(squeeze(BiexpFitArray(CurIdxs,:,3)) + CurSEADC1Array*TFactor,2);
MedianMaxADC2CIValues = median(squeeze(BiexpFitArray(CurIdxs,:,4)) + CurSEADC2Array*TFactor,2);

plot(h1,repmat(MedianMinAmp1CIValues(1),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h2,repmat(MedianMinAmp2CIValues(1),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h3,repmat(MedianMinADC1CIValues(1),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h4,repmat(MedianMinADC2CIValues(1),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h5,repmat(MedianMinAmp1CIValues(2),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h6,repmat(MedianMinAmp2CIValues(2),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h7,repmat(MedianMinADC1CIValues(2),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h8,repmat(MedianMinADC2CIValues(2),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h9,repmat(MedianMinAmp1CIValues(3),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h10,repmat(MedianMinAmp2CIValues(3),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h11,repmat(MedianMinADC1CIValues(3),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h12,repmat(MedianMinADC2CIValues(3),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');

plot(h1,repmat(MedianMaxAmp1CIValues(1),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h2,repmat(MedianMaxAmp2CIValues(1),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h3,repmat(MedianMaxADC1CIValues(1),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h4,repmat(MedianMaxADC2CIValues(1),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h5,repmat(MedianMaxAmp1CIValues(2),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h6,repmat(MedianMaxAmp2CIValues(2),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h7,repmat(MedianMaxADC1CIValues(2),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h8,repmat(MedianMaxADC2CIValues(2),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h9,repmat(MedianMaxAmp1CIValues(3),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h10,repmat(MedianMaxAmp2CIValues(3),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h11,repmat(MedianMaxADC1CIValues(3),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h12,repmat(MedianMaxADC2CIValues(3),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now repeat for a SF1 = 0.5 but different D1/D2 ratios of 8, 4, and 2
% Closest to D1/D2 = 8
find(NFTestParameterArray(:,1) > 0.497 & NFTestParameterArray(:,1) < 0.503 & 1./NFTestParameterArray(:,4) > 7.95 & 1./NFTestParameterArray(:,4) < 8.05 )
SF05_D1D28_Index = 0000;
% Closest to D1/D2 = 4
find(NFTestParameterArray(:,1) > 0.497 & NFTestParameterArray(:,1) < 0.503 & 1./NFTestParameterArray(:,4) > 3.95 & 1./NFTestParameterArray(:,4) < 4.05 )
SF05_D1D24_Index = 0000;
% Closest to D1/D2 = 2
find(NFTestParameterArray(:,1) > 0.495 & NFTestParameterArray(:,1) < 0.505 & 1./NFTestParameterArray(:,4) > 2 & 1./NFTestParameterArray(:,4) < 2.05 )
SF05_D1D22_Index = 0000;

% 3x4 Histograms
figure('Position', [50,00,1000,1000]);
h1 = subplot(4,3,1);histogram(BiexpFitArray(SF05_D1D28_Index,:,1),0:0.02:1);title('\itSF_{1}\rm\bf = 0.5, \itD_{1} / D_{2}\rm\bf = 8');ylabel('\itA_{1}\rm\bf Estimates','FontWeight','bold');ylim([0 40]);xlim([0 1]);
hold on; plot(repmat(0.5,1,2),[0 40], 'linewidth',2, 'linestyle','--'); % or 1.5*max(h.Values) if needed for adjustable
h2 = subplot(4,3,4);histogram(BiexpFitArray(SF05_D1D28_Index,:,2),0:0.02:1);ylabel('\itA_{2}\rm\bf Estimates','FontWeight','bold');ylim([0 40]);xlim([0 1]);
hold on; plot(repmat(0.5,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h3 = subplot(4,3,7);histogram(BiexpFitArray(SF05_D1D28_Index,:,3),0:0.08:4);ylabel('\itD_{1}\rm\bf Estimates','FontWeight','bold');ylim([0 40]);xlim([0 4]);
hold on; plot(repmat(1,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h4 = subplot(4,3,10);histogram(BiexpFitArray(SF05_D1D28_Index,:,4),0:0.02:1);ylabel('\itD_{2}\rm\bf Estimates','FontWeight','bold');ylim([0 40]);xlim([-0.25 0.8]);
hold on; plot(repmat(0.125,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h5 = subplot(4,3,2);histogram(BiexpFitArray(SF05_D1D24_Index,:,1),0:0.02:1);title('\itSF_{1}\rm\bf = 0.5, \itD_{1} / D_{2}\rm\bf = 4');ylim([0 40]);xlim([0 1.2]);
hold on; plot(repmat(0.5,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h6 = subplot(4,3,5);histogram(BiexpFitArray(SF05_D1D24_Index,:,2),0:0.02:1);ylim([0 40]);xlim([-0.2 1]);
hold on; plot(repmat(0.5,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h7 = subplot(4,3,8);histogram(BiexpFitArray(SF05_D1D24_Index,:,3),0:0.08:4);ylim([0 40]);xlim([0 4]);
hold on; plot(repmat(1,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h8 = subplot(4,3,11);histogram(BiexpFitArray(SF05_D1D24_Index,:,4),0:0.02:1);ylim([0 40]);xlim([-0.25 0.8]);
hold on; plot(repmat(0.25,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h9 = subplot(4,3,3);histogram(BiexpFitArray(SF05_D1D22_Index,:,1),0:0.02:1);title('\itSF_{1}\rm\bf = 0.5, \itD_{1} / D_{2}\rm\bf = 2');ylim([0 40]);xlim([0 1.2]);
hold on; plot(repmat(0.5,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h10 = subplot(4,3,6);histogram(BiexpFitArray(SF05_D1D22_Index,:,2),0:0.02:1);ylim([0 40]);xlim([-0.20 1]);
hold on; plot(repmat(0.5,1,2),[0 40], 'linewidth',2, 'linestyle','--');
h11 = subplot(4,3,9);histogram(BiexpFitArray(SF05_D1D22_Index,:,3),0:0.08:4);ylim([0 80]);xlim([0 4]);
hold on; plot(repmat(1,1,2),[0 80], 'linewidth',2, 'linestyle','--');
h12 = subplot(4,3,12);histogram(BiexpFitArray(SF05_D1D22_Index,:,4),0:0.02:1);ylim([0 80]);xlim([-0.25 0.8]);
hold on; plot(repmat(0.5,1,2),[0 80], 'linewidth',2, 'linestyle','--');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADD CONFIDENCE INTERVALS IF NEEDED
% Load the covariance matrix to extract the parameter standard errors
CurIdxs = [SF05_D1D28_Index,SF05_D1D24_Index,SF05_D1D22_Index];
CurSEAmp1Array = sqrt(CovMxArray(CurIdxs,:,1,1));
CurSEAmp2Array = sqrt(CovMxArray(CurIdxs,:,2,2));
CurSEADC1Array = sqrt(CovMxArray(CurIdxs,:,3,3));
CurSEADC2Array = sqrt(CovMxArray(CurIdxs,:,4,4));

MedianMinAmp1CIValues = median(squeeze(BiexpFitArray(CurIdxs,:,1)) - CurSEAmp1Array*TFactor,2);
MedianMinAmp2CIValues = median(squeeze(BiexpFitArray(CurIdxs,:,2)) - CurSEAmp2Array*TFactor,2);
MedianMinADC1CIValues = median(squeeze(BiexpFitArray(CurIdxs,:,3)) - CurSEADC1Array*TFactor,2);
MedianMinADC2CIValues = median(squeeze(BiexpFitArray(CurIdxs,:,4)) - CurSEADC2Array*TFactor,2);
MedianMaxAmp1CIValues = median(squeeze(BiexpFitArray(CurIdxs,:,1)) + CurSEAmp1Array*TFactor,2);
MedianMaxAmp2CIValues = median(squeeze(BiexpFitArray(CurIdxs,:,2)) + CurSEAmp2Array*TFactor,2);
MedianMaxADC1CIValues = median(squeeze(BiexpFitArray(CurIdxs,:,3)) + CurSEADC1Array*TFactor,2);
MedianMaxADC2CIValues = median(squeeze(BiexpFitArray(CurIdxs,:,4)) + CurSEADC2Array*TFactor,2);

plot(h1,repmat(MedianMinAmp1CIValues(1),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h2,repmat(MedianMinAmp2CIValues(1),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h3,repmat(MedianMinADC1CIValues(1),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h4,repmat(MedianMinADC2CIValues(1),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h5,repmat(MedianMinAmp1CIValues(2),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h6,repmat(MedianMinAmp2CIValues(2),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h7,repmat(MedianMinADC1CIValues(2),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h8,repmat(MedianMinADC2CIValues(2),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h9,repmat(MedianMinAmp1CIValues(3),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h10,repmat(MedianMinAmp2CIValues(3),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h11,repmat(MedianMinADC1CIValues(3),1,2),[0 80], 'linewidth',2, 'linestyle','--','color','b');
plot(h12,repmat(MedianMinADC2CIValues(3),1,2),[0 80], 'linewidth',2, 'linestyle','--','color','b');

plot(h1,repmat(MedianMaxAmp1CIValues(1),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h2,repmat(MedianMaxAmp2CIValues(1),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h3,repmat(MedianMaxADC1CIValues(1),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h4,repmat(MedianMaxADC2CIValues(1),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h5,repmat(MedianMaxAmp1CIValues(2),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h6,repmat(MedianMaxAmp2CIValues(2),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h7,repmat(MedianMaxADC1CIValues(2),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h8,repmat(MedianMaxADC2CIValues(2),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h9,repmat(MedianMaxAmp1CIValues(3),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h10,repmat(MedianMaxAmp2CIValues(3),1,2),[0 40], 'linewidth',2, 'linestyle','--','color','b');
plot(h11,repmat(MedianMaxADC1CIValues(3),1,2),[0 80], 'linewidth',2, 'linestyle','--','color','b');
plot(h12,repmat(MedianMaxADC2CIValues(3),1,2),[0 80], 'linewidth',2, 'linestyle','--','color','b');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
