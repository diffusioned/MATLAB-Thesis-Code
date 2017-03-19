% Chap2Display2SDMaxSignalValue.m
% MATLAB file for fitting simulated data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap2Display2SDMaxSignalValue()

load('YourPath\BiexpMaxSignalIndexArray.mat');
MedianMaxSignalArray = zeros(SF1BinDim, D1D2BinDim);

for i = 1:SF1BinDim
    for j = 1:D1D2BinDim
        CurBinIdxs = find(SF1Idxs == i & D1D2Idxs == j);
        MaxSignals = BiexpMaxSignalIndexArray(CurBinIdxs,:);
        MedianMaxSignalArray(i,j) = median(MaxSignals(:));
    end
end

PlotFontSize = 18;
figure('Position', [0,0,800,600]);
h=pcolor(SF1BinEdges(2:end), D1D2BinEdges(2:end), squeeze(MedianMaxSignalArray(:,:,1))');
% caxis([7 11]);
set(h,'edgecolor','none');
% title('Median Max Signal');
xlabel('True \itSF_{1}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
ylabel('True \itD_{1} \rm\bf/ \itD_{2}\rm\bf Value','FontWeight','bold','fontsize', PlotFontSize);
set(gca,'FontSize',PlotFontSize,'FontWeight','bold'); % Fix axis label font size
HandleC = colorbar('location', 'EastOutside','Ticks',[7,8,9,10,11]);
HandleCLabel = ylabel(HandleC,'Number Of Measurments Used');
drawnow; % MUST STICK THIS HERE TO GET THE NEXT PARTS TO WORK
LabelPos = HandleCLabel.Position;
HandleCLabel.Rotation = -90;
HandleCLabel.Position = [LabelPos(1) + 1.1, LabelPos(2), LabelPos(3)];