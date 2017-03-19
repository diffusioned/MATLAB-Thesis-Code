function Chap5DisplayBiexpFitWithNormalityBootElim()

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
% ILL-CONDITIONING
% Get the biexponential model bootstrap fits
% Reshaped earlier to get a 2D slice
load('YourPath\SliceBootFitArray.mat');
[~,~,BootDim,~] = size(SliceBiexpBootFits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now reshape the voxels to analyse the Delta AIC values and the parameter fit estimates
ReshapedBiexpBootFits = reshape(SliceBiexpBootFits,XDim*YDim,BootDim,4);
S1NormBiexpBootFits = ReshapedBiexpBootFits(HalfStrMask1 == 1,:,:);
S2NormBiexpBootFits = ReshapedBiexpBootFits(HalfStrMask2 == 1,:,:);
E1NormBiexpBootFits = ReshapedBiexpBootFits(HalfEpiMask1 == 1,:,:);
E2NormBiexpBootFits = ReshapedBiexpBootFits(HalfEpiMask2 == 1,:,:);
C1NormBiexpBootFits = ReshapedBiexpBootFits(HalfCaMask == 1,:,:);
% Save the indices from the original slice map, too
S1VoxelMaskIndices = find(HalfStrMask1 == 1);
S2VoxelMaskIndices = find(HalfStrMask2 == 1);
E1VoxelMaskIndices = find(HalfEpiMask1 == 1);
E2VoxelMaskIndices = find(HalfEpiMask2 == 1);
C1VoxelMaskIndices = find(HalfCaMask == 1);
% Also reshape fit parameters
SliceBiexpFitParams = squeeze(BiexpFitArray(:,:,17,:));
VoxelBiexpFitParams = reshape(SliceBiexpFitParams,XDim*YDim,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now identify the combined non-normal parameter estimates to correlate with ill-conditioning
AlphaSig = 0.001;
TestResults = zeros(1,4);
S1CombinedVoxelParamsTest = zeros(1, length(S1VoxelMaskIndices));
for i = 1:length(S1VoxelMaskIndices)
    % Run the lilliefors test on all parameter estimates
    for j = 1:4
        TestResults(j) = lillietest(squeeze(S1NormBiexpBootFits(i,:,j)),'Alpha',AlphaSig);
    end
    % Do a logical OR across the results and set yes for any tests that passed
    if TestResults(1) == 0 || TestResults(2) == 0 || TestResults(3) == 0 || TestResults(4) == 0
        S1CombinedVoxelParamsTest(i) = 0;
    else
        S1CombinedVoxelParamsTest(i) = 1;
    end
end

% Check count of zero.  Grab the indices that passed that test
S1NormalTestPassedIndices = S1VoxelMaskIndices(S1CombinedVoxelParamsTest == 0);
% Reshape the original parameters
S1NormalTestPassedParams = VoxelBiexpFitParams(S1NormalTestPassedIndices,:,:);
S1NormalTestPassedBootFits = S1NormBiexpBootFits(S1CombinedVoxelParamsTest == 0,:,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Also for S2
TestResults = zeros(1,4);
S2CombinedVoxelParamsTest = zeros(1, length(S2VoxelMaskIndices));
for i = 1:length(S2VoxelMaskIndices)
    % Run the lilliefors test on all parameter estimates
    for j = 1:4
        TestResults(j) = lillietest(squeeze(S2NormBiexpBootFits(i,:,j)),'Alpha',AlphaSig);
    end
    % Do a logical OR across the results and set yes for any tests that passed
    if TestResults(1) == 0 || TestResults(2) == 0 || TestResults(3) == 0 || TestResults(4) == 0
        S2CombinedVoxelParamsTest(i) = 0;
    else
        S2CombinedVoxelParamsTest(i) = 1;
    end
end

% Check count of zero.  Grab the indices that passed that test
S2NormalTestPassedIndices = S2VoxelMaskIndices(S2CombinedVoxelParamsTest == 0);
% Reshape the original parameters
S2NormalTestPassedParams = VoxelBiexpFitParams(S2NormalTestPassedIndices,:,:);
S2NormalTestPassedBootFits = S2NormBiexpBootFits(S2CombinedVoxelParamsTest == 0,:,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Also for E1
TestResults = zeros(1,4);
E1CombinedVoxelParamsTest = zeros(1, length(E1VoxelMaskIndices));
for i = 1:length(E1VoxelMaskIndices)
    % Run the lilliefors test on all parameter estimates
    for j = 1:4
        TestResults(j) = lillietest(squeeze(E1NormBiexpBootFits(i,:,j)),'Alpha',AlphaSig);
    end
    % Do a logical OR across the results and set yes for any tests that passed
    if TestResults(1) == 0 || TestResults(2) == 0 || TestResults(3) == 0 || TestResults(4) == 0
        E1CombinedVoxelParamsTest(i) = 0;
    else
        E1CombinedVoxelParamsTest(i) = 1;
    end
end

% Check count of zero.  Grab the indices that passed that test
E1NormalTestPassedIndices = E1VoxelMaskIndices(E1CombinedVoxelParamsTest == 0);
% Reshape the original parameters
E1NormalTestPassedParams = VoxelBiexpFitParams(E1NormalTestPassedIndices,:,:);
E1NormalTestPassedBootFits = E1NormBiexpBootFits(E1CombinedVoxelParamsTest == 0,:,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Also for E2
TestResults = zeros(1,4);
E2CombinedVoxelParamsTest = zeros(1, length(E2VoxelMaskIndices));
for i = 1:length(E2VoxelMaskIndices)
    % Run the lilliefors test on all parameter estimates
    for j = 1:4
        TestResults(j) = lillietest(squeeze(E2NormBiexpBootFits(i,:,j)),'Alpha',AlphaSig);
    end
    % Do a logical OR across the results and set yes for any tests that passed
    if TestResults(1) == 0 || TestResults(2) == 0 || TestResults(3) == 0 || TestResults(4) == 0
        E2CombinedVoxelParamsTest(i) = 0;
    else
        E2CombinedVoxelParamsTest(i) = 1;
    end
end

% Check count of zero.  Grab the indices that passed that test
E2NormalTestPassedIndices = E2VoxelMaskIndices(E2CombinedVoxelParamsTest == 0);
% Reshape the original parameters
E2NormalTestPassedParams = VoxelBiexpFitParams(E2NormalTestPassedIndices,:,:);
E2NormalTestPassedBootFits = E2NormBiexpBootFits(E2CombinedVoxelParamsTest == 0,:,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Also for C1
TestResults = zeros(1,4);
C1CombinedVoxelParamsTest = zeros(1, length(C1VoxelMaskIndices));
for i = 1:length(C1VoxelMaskIndices)
    % Run the lilliefors test on all parameter estimates
    for j = 1:4
        TestResults(j) = lillietest(squeeze(C1NormBiexpBootFits(i,:,j)),'Alpha',AlphaSig);
    end
    % Do a logical OR across the results and set yes for any tests that passed
    if TestResults(1) == 0 || TestResults(2) == 0 || TestResults(3) == 0 || TestResults(4) == 0
        C1CombinedVoxelParamsTest(i) = 0;
    else
        C1CombinedVoxelParamsTest(i) = 1;
    end
end

% Check count of zero.  Grab the indices that passed that test
C1NormalTestPassedIndices = C1VoxelMaskIndices(C1CombinedVoxelParamsTest == 0);
% Reshape the original parameters
C1NormalTestPassedParams = VoxelBiexpFitParams(C1NormalTestPassedIndices,:,:);
C1NormalTestPassedBootFits = C1NormBiexpBootFits(C1CombinedVoxelParamsTest == 0,:,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now display the results
% Do a 4x5 plots
figure('Position', [0,0,1200,800]);
% S1
subplot(5,4,1);histogram(squeeze(S1NormalTestPassedParams(:,1)),0:600:30000); title('\itA_{1}\rm\bf Estimates', 'fontweight','bold','fontsize',10); ylabel('S1 Voxels','FontWeight','bold','fontsize', 10);
subplot(5,4,2);histogram(squeeze(S1NormalTestPassedParams(:,2)),0:600:30000); title('\itA_{2}\rm\bf Estimates', 'fontweight','bold','fontsize',10); 
subplot(5,4,3);histogram(squeeze(1000.*S1NormalTestPassedParams(:,3)),0:0.16:8); title('\itD_{1}\rm\bf Estimates', 'fontweight','bold','fontsize',10);
subplot(5,4,4);histogram(squeeze(1000.*S1NormalTestPassedParams(:,4)),0:0.02:1); title('\itD_{2}\rm\bf Estimates', 'fontweight','bold','fontsize',10);
% S2
subplot(5,4,5);histogram(squeeze(S2NormalTestPassedParams(:,1)),0:600:30000); ylabel('S2 Voxels','FontWeight','bold','fontsize', 10);
subplot(5,4,6);histogram(squeeze(S2NormalTestPassedParams(:,2)),0:600:30000); 
subplot(5,4,7);histogram(squeeze(1000.*S2NormalTestPassedParams(:,3)),0:0.16:8); 
subplot(5,4,8);histogram(squeeze(1000.*S2NormalTestPassedParams(:,4)),0:0.02:1); 
% E1
subplot(5,4,9);histogram(squeeze(E1NormalTestPassedParams(:,1)),0:600:30000); ylabel('E1 Voxels','FontWeight','bold','fontsize', 10);
subplot(5,4,10);histogram(squeeze(E1NormalTestPassedParams(:,2)),0:600:30000);
subplot(5,4,11);histogram(squeeze(1000.*E1NormalTestPassedParams(:,3)),0:0.16:8);
subplot(5,4,12);histogram(squeeze(1000.*E1NormalTestPassedParams(:,4)),0:0.02:1);
% E2
subplot(5,4,13);histogram(squeeze(E2NormalTestPassedParams(:,1)),0:600:30000); ylabel('E2 Voxels','FontWeight','bold','fontsize', 10);
subplot(5,4,14);histogram(squeeze(E2NormalTestPassedParams(:,2)),0:600:30000); 
subplot(5,4,15);histogram(squeeze(1000.*E2NormalTestPassedParams(:,3)),0:0.16:8); 
subplot(5,4,16);histogram(squeeze(1000.*E2NormalTestPassedParams(:,4)),0:0.02:1); 
% C1
subplot(5,4,17);histogram(squeeze(C1NormalTestPassedParams(:,1)),0:600:30000); ylabel('C1 Voxels','FontWeight','bold','fontsize', 10);
subplot(5,4,18);histogram(squeeze(C1NormalTestPassedParams(:,2)),0:600:30000); 
subplot(5,4,19);histogram(squeeze(1000.*C1NormalTestPassedParams(:,3)),0:0.16:8); 
subplot(5,4,20);histogram(squeeze(1000.*C1NormalTestPassedParams(:,4)),0:0.02:1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the confidence intervals for an ROI
% Just change this line
CIPlotParamArray = C1NormalTestPassedParams;
CIPlotBootFits = C1NormalTestPassedBootFits;
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
