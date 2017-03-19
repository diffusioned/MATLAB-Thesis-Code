% Chap5DisplayBiexpFitROIBootCI.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap5DisplayBiexpFitROIBootCI()

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
% reshape the masks to pull the voxels on the 80 micron images
HalfStrMask1 = imresize(StromaMask1, 0.5);
HalfStrMask2 = imresize(StromaMask2, 0.5);
HalfEpiMask1 = imresize(EpithMask1, 0.5);
HalfEpiMask2 = imresize(EpithMask2, 0.5);
HalfCaMask = imresize(CancerMask, 0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Also bootstrap info
load('YourPath\SliceBootFitArray.mat');
[~,~,BootDim,~] = size(SliceBiexpBootFits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now reshape the voxels to analyse the Delta AIC values and the parameter fit estimates
AllVoxels = SliceBiexpBootFits(:);
S1Voxels = AllVoxels(HalfStrMask1 == 1);
S2Voxels = AllVoxels(HalfStrMask2 == 1);
E1Voxels = AllVoxels(HalfEpiMask1 == 1);
E2Voxels = AllVoxels(HalfEpiMask2 == 1);
C1Voxels = AllVoxels(HalfCaMask == 1);
% Save the indices from the original slice map, too
S1VoxelMaskIndices = find(HalfStrMask1 == 1);
S2VoxelMaskIndices = find(HalfStrMask2 == 1);
E1VoxelMaskIndices = find(HalfEpiMask1 == 1);
E2VoxelMaskIndices = find(HalfEpiMask2 == 1);
C1VoxelMaskIndices = find(HalfCaMask == 1);

% Also reshape fit parameters
SliceBiexpFitParams = squeeze(BiexpFitArray(:,:,17,:));
VoxelBiexpFitParams = reshape(SliceBiexpFitParams,XDim*YDim,4);

S1VoxelParams = VoxelBiexpFitParams(S1VoxelMaskIndices,:);
S2VoxelParams = VoxelBiexpFitParams(S2VoxelMaskIndices,:);
E1VoxelParams = VoxelBiexpFitParams(E1VoxelMaskIndices,:);
E2VoxelParams = VoxelBiexpFitParams(E2VoxelMaskIndices,:);
C1VoxelParams = VoxelBiexpFitParams(C1VoxelMaskIndices,:);

% And the bootstrap fits
ReshapedBiexpBootFits = reshape(SliceBiexpBootFits,XDim*YDim,BootDim,4);
S1BiexpBootFits = ReshapedBiexpBootFits(S1VoxelMaskIndices,:,:);
S2BiexpBootFits = ReshapedBiexpBootFits(S2VoxelMaskIndices,:,:);
E1BiexpBootFits = ReshapedBiexpBootFits(E1VoxelMaskIndices,:,:);
E2BiexpBootFits = ReshapedBiexpBootFits(E2VoxelMaskIndices,:,:);
C1BiexpBootFits = ReshapedBiexpBootFits(C1VoxelMaskIndices,:,:);

% pausehere = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the confidence intervals for an ROI
% Just change this line
CIPlotParamArray = C1VoxelParams;
CIPlotBootFits = C1BiexpBootFits;
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
%             if CurCI(2) > 50000
%                 text(50000.1, j, num2str(CurCI(2),'%1.1e'),'Color','b');
%             end
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
%             if CurCI(2) > 10
%                 text(15.1, j, num2str(CurCI(2),'%2.1f'),'Color','b');
%             end
%             if CurCI(2) > 200
%                 j
%                 tempstop=1;
%             end
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
