% Chap5DisplayHistogramsForFitTissueData.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap5DisplayHistogramsForFitTissueData()

% Load noisy data for reference
load('YourPath\YourImageData.mat');
BDim = length(BValueArray);

% Load model fits
load('YourPath\BiexpFitArray.mat');
load('YourPath\KurtFitArray.mat');
load('YourPath\MonoexpFitArray.mat');
PlotFontSize = 10;

% Set the dimension sizes.  Last dimension is the model parameters
[XDim, YDim, ZDim, ~] = size(BiexpFitArray);
TotalSignalDim = XDim*YDim*ZDim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD THE ROI DATA
load('YourPath\VoxelSelectionRegions.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reshape the masks to pull the voxels on the 80 micron images
HalfStrMask1 = imresize(StromaMask1, 0.5);
HalfStrMask2 = imresize(StromaMask2, 0.5);
HalfEpiMask1 = imresize(EpithMask1, 0.5);
HalfEpiMask2 = imresize(EpithMask2, 0.5);
HalfCaMask = imresize(CancerMask, 0.5);

% Now separate out the fitted parameter values against the masks and analyze the data
NumStromaVoxels1 = sum(HalfStrMask1(:));
NumStromaVoxels2 =  sum(HalfStrMask2(:));
NumEpithVoxels1 = sum(HalfEpiMask1(:));
NumEpithVoxels2 =  sum(HalfEpiMask2(:));
NumCancerVoxels = sum(HalfCaMask(:));

% Now reshape the voxels to analyse the parameter fit estimates
% Declare Data Arrays
S1BiexpParamEstimatesArray = zeros(XDim, YDim, 4);
S1Voxels = zeros(NumStromaVoxels1,4);
S2BiexpParamEstimatesArray = zeros(XDim, YDim, 4);
S2Voxels = zeros(NumStromaVoxels2,4);
E1BiexpParamEstimatesArray = zeros(XDim, YDim, 4);
E1Voxels = zeros(NumEpithVoxels1,4);
E2BiexpParamEstimatesArray = zeros(XDim, YDim, 4);
E2Voxels = zeros(NumEpithVoxels2,4);
C1BiexpParamEstimatesArray = zeros(XDim, YDim, 4);
C1Voxels = zeros(NumCancerVoxels,4);
for i = 1:4
    S1BiexpParamEstimatesArray(:,:,i) = squeeze(BiexpFitArray(:,:,17,i)) .* HalfStrMask1;
    CurS1Params = S1BiexpParamEstimatesArray(:,:,i);
    CurS1Params = CurS1Params(:);
    S1Voxels(:,i) = CurS1Params(HalfStrMask1 == 1);
    
    S2BiexpParamEstimatesArray(:,:,i) = squeeze(BiexpFitArray(:,:,17,i)) .* HalfStrMask2;
    CurS2Params = S2BiexpParamEstimatesArray(:,:,i);
    CurS2Params = CurS2Params(:);
    S2Voxels(:,i) = CurS2Params(HalfStrMask2 == 1);
    
    E1BiexpParamEstimatesArray(:,:,i) = squeeze(BiexpFitArray(:,:,17,i)) .* HalfEpiMask1;
    CurE1Params = E1BiexpParamEstimatesArray(:,:,i);
    CurE1Params = CurE1Params(:);
    E1Voxels(:,i) = CurE1Params(HalfEpiMask1 == 1);
    
    E2BiexpParamEstimatesArray(:,:,i) = squeeze(BiexpFitArray(:,:,17,i)) .* HalfEpiMask2;
    CurE2Params = E2BiexpParamEstimatesArray(:,:,i);
    CurE2Params = CurE2Params(:);
    E2Voxels(:,i) = CurE2Params(HalfEpiMask2 == 1);
    
    C1BiexpParamEstimatesArray(:,:,i) = squeeze(BiexpFitArray(:,:,17,i)) .* HalfCaMask;
    CurC1Params = C1BiexpParamEstimatesArray(:,:,i);
    CurC1Params = CurC1Params(:);
    C1Voxels(:,i) = CurC1Params(HalfCaMask == 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyse the voxel estimates via histogram
% Do a 4x5 plots
figure('Position', [0,0,1200,800]);
% S1
subplot(5,4,1);histogram(squeeze(S1Voxels(:,1)),0:600:30000); title('\itA_{1}\rm\bf Estimates', 'fontweight','bold','fontsize',10); ylabel('S1 Voxels','FontWeight','bold','fontsize', 10);
subplot(5,4,2);histogram(squeeze(S1Voxels(:,2)),0:600:30000); title('\itA_{2}\rm\bf Estimates', 'fontweight','bold','fontsize',10); 
subplot(5,4,3);histogram(squeeze(1000.*S1Voxels(:,3)),0:0.16:8); title('\itD_{1}\rm\bf Estimates', 'fontweight','bold','fontsize',10);
subplot(5,4,4);histogram(squeeze(1000.*S1Voxels(:,4)),0:0.02:1); title('\itD_{2}\rm\bf Estimates', 'fontweight','bold','fontsize',10);
% S2
subplot(5,4,5);histogram(squeeze(S2Voxels(:,1)),0:600:30000); ylabel('S2 Voxels','FontWeight','bold','fontsize', 10);
subplot(5,4,6);histogram(squeeze(S2Voxels(:,2)),0:600:30000); 
subplot(5,4,7);histogram(squeeze(1000.*S2Voxels(:,3)),0:0.16:8); 
subplot(5,4,8);histogram(squeeze(1000.*S2Voxels(:,4)),0:0.02:2); 
% E1
subplot(5,4,9);histogram(squeeze(E1Voxels(:,1)),0:600:30000); ylabel('E1 Voxels','FontWeight','bold','fontsize', 10);
subplot(5,4,10);histogram(squeeze(E1Voxels(:,2)),0:600:30000);
subplot(5,4,11);histogram(squeeze(1000.*E1Voxels(:,3)),0:0.16:8);
subplot(5,4,12);histogram(squeeze(1000.*E1Voxels(:,4)),0:0.02:1);
% E2
subplot(5,4,13);histogram(squeeze(E2Voxels(:,1)),0:600:30000); ylabel('E2 Voxels','FontWeight','bold','fontsize', 10);
subplot(5,4,14);histogram(squeeze(E2Voxels(:,2)),0:600:30000); 
subplot(5,4,15);histogram(squeeze(1000.*E2Voxels(:,3)),0:0.16:8); 
subplot(5,4,16);histogram(squeeze(1000.*E2Voxels(:,4)),0:0.02:1); 
% C1
subplot(5,4,17);histogram(squeeze(C1Voxels(:,1)),0:600:30000); ylabel('C1 Voxels','FontWeight','bold','fontsize', 10);
subplot(5,4,18);histogram(squeeze(C1Voxels(:,2)),0:600:30000); 
subplot(5,4,19);histogram(squeeze(1000.*C1Voxels(:,3)),0:0.16:8); 
subplot(5,4,20);histogram(squeeze(1000.*C1Voxels(:,4)),0:0.02:1); 

% Compare the mean and SD - See thesis table
mean(S1Voxels);
std(S1Voxels); 

