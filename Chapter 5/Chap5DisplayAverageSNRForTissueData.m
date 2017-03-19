% Chap5DisplayAverageSNRForTissueData.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap5DisplayAverageSNRForTissueData()
% Image Data
load('YourPath\YourImageData.mat');
[XDim, YDim, ZDim, BDim] = size(Axis1SignalArray);

% ROI Data
load('YourPath\VoxelSelectionRegions.mat');
% reshape the masks to pull the voxels on the 80 micron images
HalfStrMask1 = imresize(StromaMask1, 0.5);
HalfStrMask2 = imresize(StromaMask2, 0.5);
HalfEpiMask1 = imresize(EpithMask1, 0.5);
HalfEpiMask2 = imresize(EpithMask2, 0.5);
HalfCaMask = imresize(CancerMask, 0.5);

% Output the noise std. dev.
STDNoise
% Signal Averaged SNR Array
MeanSNRArray = zeros(XDim, YDim, ZDim); 

% SNR at lowest b value
SNRb0Array = Axis1SignalArray(:,:,:,1)./STDNoise;

% Now loop through and calculate all values for each voxel.  Flatten the 2D signal array to 1D.
for i = 1:XDim
    for j = 1:YDim
        for k = 1:ZDim
            % Divide the noise-free measured signal at each b-value by the STD of the noise to get the SNR at each point
            SNRPointArray = Axis1SignalArray(i,j,k,:)/STDNoise;
            MeanSNRArray(i,j,k) = mean(SNRPointArray);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PlotFontSize = 12;
% Colormap for b0 SNR
figure('Position', [0,0,1200,500]);subplot(2,1,1);
% h=pcolor(flipud(squeeze(MeanSNRArray(:,:,17))));set(h,'edgecolor','none');hold on;
h=imagesc(squeeze(SNRb0Array(:,:,17)));hold on;
set(gca,'xtick',[],'ytick',[],'xticklabel',{},'yticklabel',{});
% xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
% caxis([12 21]);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
% colormap(jet); % gives a good contrast for negative numbers
HandleC =  colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'SNR_{\itb\rm\bf=0}');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];

% Add contours
BndStr1 = bwboundaries(HalfStrMask1,'noholes'); boundary = BndStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2, 'LineStyle', '-')
BndStr2 = bwboundaries(HalfStrMask2,'noholes'); boundary = BndStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2, 'LineStyle', '-')
EpiStr1 = bwboundaries(HalfEpiMask1,'noholes'); boundary = EpiStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2, 'LineStyle', '-')
EpiStr2 = bwboundaries(HalfEpiMask2,'noholes'); boundary = EpiStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2, 'LineStyle', '-')
CaStr = bwboundaries(HalfCaMask,'noholes'); boundary = CaStr{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2, 'LineStyle', '-')

% Add ROI text values
text(5, -3, 'S1' ,'Color','k','FontWeight','bold');
text(15, -3, 'E1' ,'Color','k','FontWeight','bold');
text(98, -3, 'E2' ,'Color','k','FontWeight','bold');
text(134, -3, 'S2' ,'Color','k','FontWeight','bold');
text(150, -3, 'C1' ,'Color','k','FontWeight','bold');

% For the signal-averaged SNR
subplot(2,1,2);
% h=pcolor(flipud(squeeze(MeanSNRArray(:,:,17))));set(h,'edgecolor','none');hold on;
h=imagesc(squeeze(MeanSNRArray(:,:,17)));hold on;
set(gca,'xtick',[],'ytick',[],'xticklabel',{},'yticklabel',{});
% xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
% ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Ratio','FontWeight','bold','fontsize', PlotFontSize);
% caxis([12 21]);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
% colormap(jet); % gives a good contrast for negative numbers
HandleC =  colorbar('location', 'EastOutside');
HandleCLabel = ylabel(HandleC,'Signal-Averaged SNR');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];

% Add contours
BndStr1 = bwboundaries(HalfStrMask1,'noholes'); boundary = BndStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2, 'LineStyle', '-')
BndStr2 = bwboundaries(HalfStrMask2,'noholes'); boundary = BndStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2, 'LineStyle', '-')
EpiStr1 = bwboundaries(HalfEpiMask1,'noholes'); boundary = EpiStr1{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2, 'LineStyle', '-')
EpiStr2 = bwboundaries(HalfEpiMask2,'noholes'); boundary = EpiStr2{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2, 'LineStyle', '-')
CaStr = bwboundaries(HalfCaMask,'noholes'); boundary = CaStr{1};
plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2, 'LineStyle', '-')

% b0 SNR values for table
CurS1b0SNR = squeeze(SNRb0Array(:,:,17)) .* HalfStrMask1;
CurS1b0SNR = CurS1b0SNR(:);
S1b0SNR = CurS1b0SNR(HalfStrMask1 == 1);
mean(S1b0SNR)

CurS2b0SNR = squeeze(SNRb0Array(:,:,17)) .* HalfStrMask2;
CurS2b0SNR = CurS2b0SNR(:);
S2b0SNR = CurS2b0SNR(HalfStrMask2 == 1);
mean(S2b0SNR)

CurE1b0SNR = squeeze(SNRb0Array(:,:,17)) .* HalfEpiMask1;
CurE1b0SNR = CurE1b0SNR(:);
E1b0SNR = CurE1b0SNR(HalfEpiMask1 == 1);
mean(E1b0SNR)

CurE2b0SNR = squeeze(SNRb0Array(:,:,17)) .* HalfEpiMask2;
CurE2b0SNR = CurE2b0SNR(:);
E2b0SNR = CurE2b0SNR(HalfEpiMask2 == 1);
mean(E2b0SNR)

CurC1b0SNR = squeeze(SNRb0Array(:,:,17)) .* HalfCaMask;
CurC1b0SNR = CurC1b0SNR(:);
C1b0SNR = CurC1b0SNR(HalfCaMask == 1);
mean(C1b0SNR)

% Average values for table
CurS1SigAvSNR = squeeze(MeanSNRArray(:,:,17)) .* HalfStrMask1;
CurS1SigAvSNR = CurS1SigAvSNR(:);
S1SigAvSNR = CurS1SigAvSNR(HalfStrMask1 == 1);
mean(S1SigAvSNR)

CurS2SigAvSNR = squeeze(MeanSNRArray(:,:,17)) .* HalfStrMask2;
CurS2SigAvSNR = CurS2SigAvSNR(:);
S2SigAvSNR = CurS2SigAvSNR(HalfStrMask2 == 1);
mean(S2SigAvSNR)

CurE1SigAvSNR = squeeze(MeanSNRArray(:,:,17)) .* HalfEpiMask1;
CurE1SigAvSNR = CurE1SigAvSNR(:);
E1SigAvSNR = CurE1SigAvSNR(HalfEpiMask1 == 1);
mean(E1SigAvSNR)

CurE2SigAvSNR = squeeze(MeanSNRArray(:,:,17)) .* HalfEpiMask2;
CurE2SigAvSNR = CurE2SigAvSNR(:);
E2SigAvSNR = CurE2SigAvSNR(HalfEpiMask2 == 1);
mean(E2SigAvSNR)

CurC1SigAvSNR = squeeze(MeanSNRArray(:,:,17)) .* HalfCaMask;
CurC1SigAvSNR = CurC1SigAvSNR(:);
C1SigAvSNR = CurC1SigAvSNR(HalfCaMask == 1);
mean(C1SigAvSNR)