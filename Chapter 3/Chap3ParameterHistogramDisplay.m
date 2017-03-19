% Chap3ParameterHistogramDisplay.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap3ParameterHistogramDisplay()

% This function is going to determine the test points to display the parameter histogram fits

% First test set - D1/D2 = 15, SF1 = 0.05, 0.1, 0.25, 0.5
load('YourPath\NoiseFreeParameters.mat');

% Note, the values below were unique to my test set.  Use the commented find statements to find similar values in yours
% Closest to SF = 0.5
% find(NFTestParameterArray(:,1) > 0.501 & NFTestParameterArray(:,1) < 0.502 & 1./NFTestParameterArray(:,4) > 14.9 & 1./NFTestParameterArray(:,4) < 15.1 )
SF05_D1D215_Index = 31849;
% Closest to SF = 0.25
% find(NFTestParameterArray(:,1) > 0.249 & NFTestParameterArray(:,1) < 0.251 & 1./NFTestParameterArray(:,4) > 14.9 & 1./NFTestParameterArray(:,4) < 15.1 )
SF025_D1D215_Index = 16282;
% Closest to SF = 0.1
% find(NFTestParameterArray(:,1) > 0.095 & NFTestParameterArray(:,1) < 0.105 & 1./NFTestParameterArray(:,4) > 14.95 & 1./NFTestParameterArray(:,4) < 15.05 )
SF01_D1D215_Index = 5522;
% Closest to SF = 0.05
% find(NFTestParameterArray(:,1) > 0.045 & NFTestParameterArray(:,1) < 0.055 & 1./NFTestParameterArray(:,4) > 14.95 & 1./NFTestParameterArray(:,4) < 15.05 )
SF005_D1D215_Index = 28156;
% Closest to SF = 0.95
% find(NFTestParameterArray(:,1) > 0.945 & NFTestParameterArray(:,1) < 0.955 & 1./NFTestParameterArray(:,4) > 14.95 & 1./NFTestParameterArray(:,4) < 15.05 )
SF095_D1D215_Index = 5811;

% Load the corresponding parameter fits
load('YourPath\KurtFitArray.mat');

CurIdxs = [31849,28156,5811];
CurSEAmpArray = sqrt(CovMxArray(CurIdxs,:,1,1));
CurSEDappArray = sqrt(CovMxArray(CurIdxs,:,2,2));
CurSEKappArray = sqrt(CovMxArray(CurIdxs,:,3,3));

% 3x3 Histograms
figure('Position', [0,0,1000,1000]);
h1 = subplot(3,3,1);histogram(KurtFitArray(SF05_D1D215_Index,:,1),0.85:0.005:1.05);title('\itSF_{1}\rm\bf = 0.5, \itD_{1} / D_{2}\rm\bf = 15');ylabel('\itS_{0}\rm\bf Estimates','FontWeight','bold');xlim([0.85 1.05]);%ylim([0 30]);
h2 = subplot(3,3,4);histogram(KurtFitArray(SF05_D1D215_Index,:,2),0.:0.01:0.6);ylabel('\itD_{app}\rm\bf Estimates','FontWeight','bold');xlim([0. 0.6]);ylim([0 40]);
h3 = subplot(3,3,7);histogram(KurtFitArray(SF05_D1D215_Index,:,3),-1:0.05:1.2);ylabel('\itK_{app}\rm\bf Estimates','FontWeight','bold');xlim([-1 1.2]);ylim([0 100]);
h4 = subplot(3,3,2);histogram(KurtFitArray(SF095_D1D215_Index,:,1),0.9:0.005:1.1);title('\itSF_{1}\rm\bf = 0.95, \itD_{1} / D_{2}\rm\bf = 15');xlim([0.9 1.1]);ylim([0 30]);
h5 = subplot(3,3,5);histogram(KurtFitArray(SF095_D1D215_Index,:,2),0.6:0.01:1.2);xlim([0.6 1.2]);ylim([0 20]);
h6 = subplot(3,3,8);histogram(KurtFitArray(SF095_D1D215_Index,:,3),-1:0.05:1.2);xlim([-1 1.2]);ylim([0 160]);
h7 = subplot(3,3,3);histogram(KurtFitArray(SF005_D1D215_Index,:,1),0.9:0.005:1.1);title('\itSF_{1}\rm\bf = 0.05, \itD_{1} / D_{2}\rm\bf = 15');xlim([0.9 1.1]);ylim([0 30]);
h8 = subplot(3,3,6);histogram(KurtFitArray(SF005_D1D215_Index,:,2),0:0.01:0.6);xlim([0. 0.6]);ylim([0 80]);
h9 = subplot(3,3,9);histogram(KurtFitArray(SF005_D1D215_Index,:,3),-1:0.05:1.2);xlim([-1 1.2]);ylim([0 20]);
