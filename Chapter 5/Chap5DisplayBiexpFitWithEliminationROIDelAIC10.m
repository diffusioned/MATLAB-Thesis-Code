% Chap5DisplayBiexpFitWithEliminationROIDelAIC10.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap5DisplayBiexpFitWithEliminationROIDelAIC10()

% Load noisy data for reference
load('YourPath\E24Axis1SignalArray.mat');
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

% reshape the masks to pull the voxels on the 80 micron images
HalfStrMask1 = imresize(StromaMask1, 0.5);
HalfStrMask2 = imresize(StromaMask2, 0.5);
HalfEpiMask1 = imresize(EpithMask1, 0.5);
HalfEpiMask2 = imresize(EpithMask2, 0.5);
HalfCaMask = imresize(CancerMask, 0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL SELECTION AND ILL-CONDITIONING
% Grab the RSS values and create a 2D slice
load('C:\Users\ccha4217\Documents\Aux Thesis Fits\E24 Tissue Testing\Diagnostics\BiexpRSSArray.mat');
load('C:\Users\ccha4217\Documents\Aux Thesis Fits\E24 Tissue Testing\Diagnostics\MonoexpRSSArray.mat');

SliceBiexpRSS = squeeze(BiexpRSSArray(:,:,17));
MonoexpParam = 2+1; KurtParam = 3+1;  BiexpParam = 4+1; % The +1 is to include the variance needed with RSS fits
BiexpAICArray = BDim.*log(BiexpRSSArray/BDim) + 2*(BiexpParam);
SliceBiexpAIC = squeeze(BiexpAICArray(:,:,17));
MonoexpAICArray = BDim.*log(MonoexpRSSArray/BDim) + 2*(MonoexpParam);
SliceMonoexpAIC = squeeze(MonoexpAICArray(:,:,17));
SliceDeltaAICBvM = SliceBiexpAIC-SliceMonoexpAIC;
% Load the boot fits
load('C:\Users\ccha4217\Documents\Aux Thesis Fits\E24 Tissue Testing\SliceBootFitArray.mat');
[~,~,BootDim,~] = size(SliceBiexpBootFits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now reshape the voxels to analyse the Delta AIC values and the parameter fit estimates
AllVoxelsDeltaAICBvM = SliceDeltaAICBvM(:);
S1AICBvMVoxels = AllVoxelsDeltaAICBvM(HalfStrMask1 == 1);
S2AICBvMVoxels = AllVoxelsDeltaAICBvM(HalfStrMask2 == 1);
E1AICBvMVoxels = AllVoxelsDeltaAICBvM(HalfEpiMask1 == 1);
E2AICBvMVoxels = AllVoxelsDeltaAICBvM(HalfEpiMask2 == 1);
C1AICBvMVoxels = AllVoxelsDeltaAICBvM(HalfCaMask == 1);
% Save the indices from the original slice map, too
S1VoxelMaskIndices = find(HalfStrMask1 == 1);
S2VoxelMaskIndices = find(HalfStrMask2 == 1);
E1VoxelMaskIndices = find(HalfEpiMask1 == 1);
E2VoxelMaskIndices = find(HalfEpiMask2 == 1);
C1VoxelMaskIndices = find(HalfCaMask == 1);

% Histograms of Delta AIC values for each region
figure('Position', [0,0,1200,800]);
subplot(3,2,1);histogram(S1AICBvMVoxels,50);
title('S1 Voxels', 'fontweight','bold','fontsize',10); xlabel('\DeltaAIC (Biexp - Monoexp)','FontWeight','bold','fontsize', 10);ylabel('Count','FontWeight','bold','fontsize', 10);
subplot(3,2,2);histogram(S2AICBvMVoxels,50);
title('S2 Voxels', 'fontweight','bold','fontsize',10); xlabel('\DeltaAIC (Biexp - Monoexp)','FontWeight','bold','fontsize', 10);ylabel('Count','FontWeight','bold','fontsize', 10);
subplot(3,2,3);histogram(E1AICBvMVoxels,50);
title('E1 Voxels', 'fontweight','bold','fontsize',10); xlabel('\DeltaAIC (Biexp - Monoexp)','FontWeight','bold','fontsize', 10);ylabel('Count','FontWeight','bold','fontsize', 10);
subplot(3,2,4);histogram(E2AICBvMVoxels,50);
title('E2 Voxels', 'fontweight','bold','fontsize',10); xlabel('\DeltaAIC (Biexp - Monoexp)','FontWeight','bold','fontsize', 10);ylabel('Count','FontWeight','bold','fontsize', 10);
subplot(3,2,5);histogram(C1AICBvMVoxels,50);
title('C1 Voxels', 'fontweight','bold','fontsize',10); xlabel('\DeltaAIC (Biexp - Monoexp)','FontWeight','bold','fontsize', 10);ylabel('Count','FontWeight','bold','fontsize', 10);

% Now remove all voxels with Delta AIC > -10
MinDeltaS1AICBvMVoxels = S1AICBvMVoxels(S1AICBvMVoxels <= -10);
MinDeltaS2AICBvMVoxels = S2AICBvMVoxels(S2AICBvMVoxels <= -10);
MinDeltaE1AICBvMVoxels = E1AICBvMVoxels(E1AICBvMVoxels <= -10);
MinDeltaE2AICBvMVoxels = E2AICBvMVoxels(E2AICBvMVoxels <= -10);
MinDeltaC1AICBvMVoxels = C1AICBvMVoxels(C1AICBvMVoxels <= -10);
% And the same with the list of indices
S1VoxelMaskIndices = S1VoxelMaskIndices(S1AICBvMVoxels <= -10);
S2VoxelMaskIndices = S2VoxelMaskIndices(S2AICBvMVoxels <= -10);
E1VoxelMaskIndices = E1VoxelMaskIndices(E1AICBvMVoxels <= -10);
E2VoxelMaskIndices = E2VoxelMaskIndices(E2AICBvMVoxels <= -10);
C1VoxelMaskIndices = C1VoxelMaskIndices(C1AICBvMVoxels <= -10);

% Also reshape fit parameters
SliceBiexpFitParams = squeeze(BiexpFitArray(:,:,17,:));
VoxelBiexpFitParams = reshape(SliceBiexpFitParams,XDim*YDim,4);
% Data Arrays for voxel parameters
S1AICBvMVoxelParams = VoxelBiexpFitParams(S1VoxelMaskIndices,:);
S2AICBvMVoxelParams = VoxelBiexpFitParams(S2VoxelMaskIndices,:);
E1AICBvMVoxelParams = VoxelBiexpFitParams(E1VoxelMaskIndices,:);
E2AICBvMVoxelParams = VoxelBiexpFitParams(E2VoxelMaskIndices,:);
C1AICBvMVoxelParams = VoxelBiexpFitParams(C1VoxelMaskIndices,:);

% Reshape boot fits, too
ReshapedBiexpBootFits = reshape(SliceBiexpBootFits,XDim*YDim,BootDim,4);
S1AICBvMBiexpBootFits = ReshapedBiexpBootFits(S1VoxelMaskIndices,:,:);
S2AICBvMBiexpBootFits = ReshapedBiexpBootFits(S2VoxelMaskIndices,:,:);
E1AICBvMBiexpBootFits = ReshapedBiexpBootFits(E1VoxelMaskIndices,:,:);
E2AICBvMBiexpBootFits = ReshapedBiexpBootFits(E2VoxelMaskIndices,:,:);
C1AICBvMBiexpBootFits = ReshapedBiexpBootFits(C1VoxelMaskIndices,:,:);

pausehere = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyse the voxel estimates via histogram
% Do a 4x5 plot
figure('Position', [0,0,1200,800]);
% S1
subplot(5,4,1);histogram(squeeze(S1AICBvMVoxelParams(:,1)),0:600:30000); title('\itA_{1}\rm\bf Estimates', 'fontweight','bold','fontsize',10); ylabel('S1 Voxels','FontWeight','bold','fontsize', 10);
subplot(5,4,2);histogram(squeeze(S1AICBvMVoxelParams(:,2)),0:600:30000); title('\itA_{2}\rm\bf Estimates', 'fontweight','bold','fontsize',10); 
subplot(5,4,3);histogram(squeeze(1000.*S1AICBvMVoxelParams(:,3)),0:0.16:8); title('\itD_{1}\rm\bf Estimates', 'fontweight','bold','fontsize',10);
subplot(5,4,4);histogram(squeeze(1000.*S1AICBvMVoxelParams(:,4)),0:0.02:1); title('\itD_{2}\rm\bf Estimates', 'fontweight','bold','fontsize',10);
% S2
subplot(5,4,5);histogram(squeeze(S2AICBvMVoxelParams(:,1)),0:600:30000); ylabel('S2 Voxels','FontWeight','bold','fontsize', 10);
subplot(5,4,6);histogram(squeeze(S2AICBvMVoxelParams(:,2)),0:600:30000); 
subplot(5,4,7);histogram(squeeze(1000.*S2AICBvMVoxelParams(:,3)),0:0.16:8); 
subplot(5,4,8);histogram(squeeze(1000.*S2AICBvMVoxelParams(:,4)),0:0.02:1); 
% E1
subplot(5,4,9);histogram(squeeze(E1AICBvMVoxelParams(:,1)),0:600:30000); ylabel('E1 Voxels','FontWeight','bold','fontsize', 10);
subplot(5,4,10);histogram(squeeze(E1AICBvMVoxelParams(:,2)),0:600:30000);
subplot(5,4,11);histogram(squeeze(1000.*E1AICBvMVoxelParams(:,3)),0:0.16:8);
subplot(5,4,12);histogram(squeeze(1000.*E1AICBvMVoxelParams(:,4)),0:0.02:1);
% E2
subplot(5,4,13);histogram(squeeze(E2AICBvMVoxelParams(:,1)),0:600:30000); ylabel('E2 Voxels','FontWeight','bold','fontsize', 10);
subplot(5,4,14);histogram(squeeze(E2AICBvMVoxelParams(:,2)),0:600:30000); 
subplot(5,4,15);histogram(squeeze(1000.*E2AICBvMVoxelParams(:,3)),0:0.16:8); 
subplot(5,4,16);histogram(squeeze(1000.*E2AICBvMVoxelParams(:,4)),0:0.02:1); 
% C1
subplot(5,4,17);histogram(squeeze(C1AICBvMVoxelParams(:,1)),0:600:30000); ylabel('C1 Voxels','FontWeight','bold','fontsize', 10);
subplot(5,4,18);histogram(squeeze(C1AICBvMVoxelParams(:,2)),0:600:30000); 
subplot(5,4,19);histogram(squeeze(1000.*C1AICBvMVoxelParams(:,3)),0:0.16:8); 
subplot(5,4,20);histogram(squeeze(1000.*C1AICBvMVoxelParams(:,4)),0:0.02:1); 

% Compare the mean and SD
DisplayAICBvMVoxelParams = C1AICBvMVoxelParams;
[num2str(mean(DisplayAICBvMVoxelParams(:,1))), ' + ', num2str(std(DisplayAICBvMVoxelParams(:,1)))]
[num2str(mean(DisplayAICBvMVoxelParams(:,2))), ' + ', num2str(std(DisplayAICBvMVoxelParams(:,2)))]
[num2str(mean(DisplayAICBvMVoxelParams(:,1)./(DisplayAICBvMVoxelParams(:,1)+DisplayAICBvMVoxelParams(:,2)))), ' + ', num2str(std(DisplayAICBvMVoxelParams(:,1)./(DisplayAICBvMVoxelParams(:,1)+DisplayAICBvMVoxelParams(:,2))))]
[num2str(mean(DisplayAICBvMVoxelParams(:,3))*1000), ' + ', num2str(std(DisplayAICBvMVoxelParams(:,3))*1000)]
[num2str(mean(DisplayAICBvMVoxelParams(:,4))*1000), ' + ', num2str(std(DisplayAICBvMVoxelParams(:,4))*1000)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the confidence intervals for an ROI
% Just change this line
CIPlotParamArray = C1AICBvMVoxelParams;
CIPlotBootFits = C1AICBvMBiexpBootFits;
NumCIVoxels = length(CIPlotParamArray(:,1));
figure('Position', [0,0,800,800]);
for i = 1:4
    subplot(4,1,i);
    for j = 1:length(CIPlotParamArray(:,1))
        if i == 1
            CurCI = prctile(squeeze(CIPlotBootFits(j,:,i)),[2.5 97.5]);
            plot([j,j], CurCI, 'b');hold on;
            set(gca, 'Ticklength', [0 0],'XTick',[]); % Inside tick marks are annoying here
            scatter(j,CIPlotParamArray(j,i),10,'r','filled','d');
            ylim([0 50000]);   xlim([0 NumCIVoxels+1]);
        elseif i == 2
            CurCI = prctile(squeeze(CIPlotBootFits(j,:,i)),[2.5 97.5]);
            plot([j,j], CurCI, 'b');hold on;
            set(gca, 'Ticklength', [0 0],'XTick',[]);
            scatter(j,CIPlotParamArray(j,i),10,'r','filled','d');
            ylim([0 30000]);   xlim([0 NumCIVoxels+1]);
        elseif i == 3
            CurCI = prctile(1000.*squeeze(CIPlotBootFits(j,:,i)),[2.5 97.5]);
            plot([j,j], CurCI, 'b');hold on;
            set(gca, 'Ticklength', [0 0],'XTick',[]);
            scatter(j,1000.*CIPlotParamArray(j,i),10,'r','filled','d');
            ylim([0 20]);   xlim([0 NumCIVoxels+1]);
        elseif i == 4
            CurCI = prctile(1000.*squeeze(CIPlotBootFits(j,:,i)),[2.5 97.5]);
            plot([j,j], CurCI, 'b');hold on;
            set(gca, 'Ticklength', [0 0],'XTick',[]);
            scatter(j,1000.*CIPlotParamArray(j,i),10,'r','filled','d');
            ylim([0 1]);   xlim([0 NumCIVoxels+1]);
        end
    end
    if i == 1
        ylabel('\itA\rm\bf_{1}','fontweight','bold');
    elseif i == 2
        ylabel('\itA\rm\bf_{2}','fontweight','bold');
    elseif i == 3
        ylabel('\itD\rm\bf_{1}','fontweight','bold');
    elseif i == 4
        ylabel('\itD\rm\bf_{2}','fontweight','bold');
    end
end
