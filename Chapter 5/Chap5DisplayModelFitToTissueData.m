% Chap5DisplayModelFitToTissueData.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap5DisplayModelFitToTissueData()

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

% This ROI data below was done in previous code and saved as a .mat file.  This can be done with this code snippet:
% imagesc(im) % show the image
% h = imfreehand; % trace the ROI
% mask = createMask(h); % create a mask from it

% ROI Data
load('YourPath\VoxelSelectionRegions.mat');
% reshape the masks to pull the voxels on the 80 micron images
HalfStrMask1 = imresize(StromaMask1, 0.5);
HalfStrMask2 = imresize(StromaMask2, 0.5);
HalfEpiMask1 = imresize(EpithMask1, 0.5);
HalfEpiMask2 = imresize(EpithMask2, 0.5);
HalfCaMask = imresize(CancerMask, 0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display amplitude parameter maps.  Note, the data to be displayed needs to be 2D, so reshape and squeeze as needed.
PlotFontSize = 12;
figure('Position', [0,0,800,800]);
subplot(5,1,1);h=imagesc(squeeze(MonoexpFitArray(:,:,17,1)));hold on;
% title('Monoexp \itA_{0}\rm\bf Fits');
set(gca,'xtick',[],'ytick',[],'xticklabel',{},'yticklabel',{});
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC = colorbar('location', 'EastOutside');caxis([0 30000]);
HandleCLabel = ylabel(HandleC,'\itS\rm\bf_{0} (mono)');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];
% Add contours
BndStr1 = bwboundaries(HalfStrMask1,'noholes'); boundary = BndStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
BndStr2 = bwboundaries(HalfStrMask2,'noholes'); boundary = BndStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr1 = bwboundaries(HalfEpiMask1,'noholes'); boundary = EpiStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr2 = bwboundaries(HalfEpiMask2,'noholes'); boundary = EpiStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
CaStr = bwboundaries(HalfCaMask,'noholes'); boundary = CaStr{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')

% Add ROI text values
text(5, -5, 'S1' ,'Color','k','FontWeight','bold');
text(15, -5, 'E1' ,'Color','k','FontWeight','bold');
text(98, -5, 'E2' ,'Color','k','FontWeight','bold');
text(134, -5, 'S2' ,'Color','k','FontWeight','bold');
text(150, -5, 'C1' ,'Color','k','FontWeight','bold');

% Kurtosis
subplot(5,1,2);h=imagesc(squeeze(KurtFitArray(:,:,17,1)));hold on;
% title('Kurtosis \itA_{0}\rm\bf Fits');
set(gca,'xtick',[],'ytick',[],'xticklabel',{},'yticklabel',{});
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC = colorbar('location', 'EastOutside');caxis([0 30000]);
HandleCLabel = ylabel(HandleC,'\itS\rm\bf_{0} (kurt)');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];
% Add contours
BndStr1 = bwboundaries(HalfStrMask1,'noholes'); boundary = BndStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
BndStr2 = bwboundaries(HalfStrMask2,'noholes'); boundary = BndStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr1 = bwboundaries(HalfEpiMask1,'noholes'); boundary = EpiStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr2 = bwboundaries(HalfEpiMask2,'noholes'); boundary = EpiStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
CaStr = bwboundaries(HalfCaMask,'noholes'); boundary = CaStr{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')

% Biexp
subplot(5,1,3);h=imagesc(squeeze(BiexpFitArray(:,:,17,1)));hold on;
% title('Biexp \itA_{1}\rm\bf Fits');
set(gca,'xtick',[],'ytick',[],'xticklabel',{},'yticklabel',{});
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC = colorbar('location', 'EastOutside');caxis([0 30000]);
HandleCLabel = ylabel(HandleC,'\itA\rm\bf_{1} (biexp)');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];
% Add contours
BndStr1 = bwboundaries(HalfStrMask1,'noholes'); boundary = BndStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
BndStr2 = bwboundaries(HalfStrMask2,'noholes'); boundary = BndStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr1 = bwboundaries(HalfEpiMask1,'noholes'); boundary = EpiStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr2 = bwboundaries(HalfEpiMask2,'noholes'); boundary = EpiStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
CaStr = bwboundaries(HalfCaMask,'noholes'); boundary = CaStr{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')

subplot(5,1,4);h=imagesc(squeeze(BiexpFitArray(:,:,17,2)));hold on;
% title('Biexp \itA_{2}\rm\bf Fits');
set(gca,'xtick',[],'ytick',[],'xticklabel',{},'yticklabel',{});
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC = colorbar('location', 'EastOutside');caxis([0 30000]);
HandleCLabel = ylabel(HandleC,'\itA\rm\bf_{2} (biexp)');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];
% Add contours
BndStr1 = bwboundaries(HalfStrMask1,'noholes'); boundary = BndStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
BndStr2 = bwboundaries(HalfStrMask2,'noholes'); boundary = BndStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr1 = bwboundaries(HalfEpiMask1,'noholes'); boundary = EpiStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr2 = bwboundaries(HalfEpiMask2,'noholes'); boundary = EpiStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
CaStr = bwboundaries(HalfCaMask,'noholes'); boundary = CaStr{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')

% SF1
subplot(5,1,5);h=imagesc(squeeze(BiexpFitArray(:,:,17,1)./(BiexpFitArray(:,:,17,1)+BiexpFitArray(:,:,17,2))));hold on;
% title('Biexp \itA_{2}\rm\bf Fits');
set(gca,'xtick',[],'ytick',[],'xticklabel',{},'yticklabel',{});
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC = colorbar('location', 'EastOutside');caxis([0 1]);
HandleCLabel = ylabel(HandleC,'\itSF\rm\bf_{1} (biexp)');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 0.6, LabelPos(2), LabelPos(3)];
% Add contours
BndStr1 = bwboundaries(HalfStrMask1,'noholes'); boundary = BndStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
BndStr2 = bwboundaries(HalfStrMask2,'noholes'); boundary = BndStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr1 = bwboundaries(HalfEpiMask1,'noholes'); boundary = EpiStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr2 = bwboundaries(HalfEpiMask2,'noholes'); boundary = EpiStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
CaStr = bwboundaries(HalfCaMask,'noholes'); boundary = CaStr{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')

pausehere = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display decay parameter maps
PlotFontSize = 10;
figure('Position', [0,0,800,800]);
subplot(5,1,1);h=imagesc(squeeze(1000*MonoexpFitArray(:,:,17,2)));hold on;
set(gca,'xtick',[],'ytick',[],'xticklabel',{},'yticklabel',{});
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC = colorbar('location', 'EastOutside');caxis([0 2.5]);
HandleCLabel = ylabel(HandleC,'\itADC\rm\bf (mono)');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];
% Add contours
BndStr1 = bwboundaries(HalfStrMask1,'noholes'); boundary = BndStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
BndStr2 = bwboundaries(HalfStrMask2,'noholes'); boundary = BndStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr1 = bwboundaries(HalfEpiMask1,'noholes'); boundary = EpiStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr2 = bwboundaries(HalfEpiMask2,'noholes'); boundary = EpiStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
CaStr = bwboundaries(HalfCaMask,'noholes'); boundary = CaStr{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')

% Add ROI text values
text(5, -5, 'S1' ,'Color','k','FontWeight','bold');
text(15, -5, 'E1' ,'Color','k','FontWeight','bold');
text(98, -5, 'E2' ,'Color','k','FontWeight','bold');
text(134, -5, 'S2' ,'Color','k','FontWeight','bold');
text(150, -5, 'C1' ,'Color','k','FontWeight','bold');

% Kurtosis
subplot(5,1,2);h=imagesc(squeeze(1000*KurtFitArray(:,:,17,2)));hold on;
set(gca,'xtick',[],'ytick',[],'xticklabel',{},'yticklabel',{});
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC = colorbar('location', 'EastOutside');caxis([0 2.5]);
HandleCLabel = ylabel(HandleC,'\itD\rm\bf_{app} (kurt)');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];
% Add contours
BndStr1 = bwboundaries(HalfStrMask1,'noholes'); boundary = BndStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
BndStr2 = bwboundaries(HalfStrMask2,'noholes'); boundary = BndStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr1 = bwboundaries(HalfEpiMask1,'noholes'); boundary = EpiStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr2 = bwboundaries(HalfEpiMask2,'noholes'); boundary = EpiStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
CaStr = bwboundaries(HalfCaMask,'noholes'); boundary = CaStr{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')

subplot(5,1,3);h=imagesc(squeeze(KurtFitArray(:,:,17,3)));hold on;
% title('Biexp \itA_{1}\rm\bf Fits');
set(gca,'xtick',[],'ytick',[],'xticklabel',{},'yticklabel',{});
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC = colorbar('location', 'EastOutside');caxis([-1 2]);
HandleCLabel = ylabel(HandleC,'\itK\rm\bf_{app} (kurt)');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];
% Add contours
BndStr1 = bwboundaries(HalfStrMask1,'noholes'); boundary = BndStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
BndStr2 = bwboundaries(HalfStrMask2,'noholes'); boundary = BndStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr1 = bwboundaries(HalfEpiMask1,'noholes'); boundary = EpiStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr2 = bwboundaries(HalfEpiMask2,'noholes'); boundary = EpiStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
CaStr = bwboundaries(HalfCaMask,'noholes'); boundary = CaStr{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')

% Biexp
subplot(5,1,4);h=imagesc(squeeze(1000*BiexpFitArray(:,:,17,3)));hold on;
% title('Biexp \itA_{2}\rm\bf Fits');
set(gca,'xtick',[],'ytick',[],'xticklabel',{},'yticklabel',{});
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC = colorbar('location', 'EastOutside');caxis([0 4]);
HandleCLabel = ylabel(HandleC,'\itD\rm\bf_{1} (biexp)');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];
% Add contours
BndStr1 = bwboundaries(HalfStrMask1,'noholes'); boundary = BndStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
BndStr2 = bwboundaries(HalfStrMask2,'noholes'); boundary = BndStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr1 = bwboundaries(HalfEpiMask1,'noholes'); boundary = EpiStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr2 = bwboundaries(HalfEpiMask2,'noholes'); boundary = EpiStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
CaStr = bwboundaries(HalfCaMask,'noholes'); boundary = CaStr{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')

subplot(5,1,5);h=imagesc(squeeze(1000*BiexpFitArray(:,:,17,4)));hold on;
% title('Biexp \itA_{2}\rm\bf Fits');
set(gca,'xtick',[],'ytick',[],'xticklabel',{},'yticklabel',{});
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC = colorbar('location', 'EastOutside');caxis([0 2]);
HandleCLabel = ylabel(HandleC,'\itD\rm\bf_{2} (biexp)');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];
% Add contours
BndStr1 = bwboundaries(HalfStrMask1,'noholes'); boundary = BndStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
BndStr2 = bwboundaries(HalfStrMask2,'noholes'); boundary = BndStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr1 = bwboundaries(HalfEpiMask1,'noholes'); boundary = EpiStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
EpiStr2 = bwboundaries(HalfEpiMask2,'noholes'); boundary = EpiStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
CaStr = bwboundaries(HalfCaMask,'noholes'); boundary = CaStr{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 1.5, 'LineStyle', '-')
