% Chap2Plot95PctGaussianConfidenceInterval.m
% MATLAB file for fitting simulated data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap2Plot95PctGaussianConfidenceInterval()

% This function is for creating a log-lin plot of the range of signals that will be teste
x = -4:.1:4;
norm = normpdf(x,0,1);

% Plot the pdf.
figure('Position', [00,0,1000,600]);
plot(x,norm,'LineWidth',4,'Color','b');hold on;
PlotFontSize = 14;
set(gca,'FontSize',PlotFontSize,'FontWeight','bold');
ylabel('Probability Density','FontWeight','bold','fontsize', PlotFontSize);
xlabel('x','FontWeight','bold','fontsize', PlotFontSize);

% Plot the 95% confidence intervals, too
hArea = area(x(21:61),norm(21:61));
% these statements don't work:
% hFace = get(hArea,'children');
% set(hFace,'FaceAlpha',0.5);
plot(repmat(-2,1,51),0:0.008:0.4,'r--','LineWidth',3);
plot(repmat(2,1,51),0:0.008:0.4,'r--','LineWidth',3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add some text explanation for all added labels/points
% text(18,0.035,'SNR 25 Noise Average','Color','k','FontSize',24 , 'FontWeight', 'bold');
stophere = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%