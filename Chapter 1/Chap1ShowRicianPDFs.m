% Chap1ShowRicianPDFs.m
% MATLAB file for creating Figure 2 in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
SNRToTest = [0.0, 2.0, 5.0];

NumberOfPDFs = length(SNRToTest);
XIncrement = 0.1;
XValue = 0:XIncrement:10;
XDim = length(XValue);
STDNoise = 1.0;

PDFArray = zeros(NumberOfPDFs,XDim);

% Create the Rician PDFs
for i = 1:NumberOfPDFs
    PDFArray(i,:) =  pdf('Rician', XValue, SNRToTest(i), STDNoise);
end
   
% Display results
figure('position',[50 50 800 600]);plot(XValue,PDFArray(1,:),'Color','r','LineWidth',2);hold all;xlim([0 10]);
plot(XValue,PDFArray(2,:),'Color','g','LineWidth',2);
plot(XValue,PDFArray(3,:),'Color','b','LineWidth',2);
set(gca,'FontSize',18)
xlabel('M/\sigma','fontsize',18, 'fontweight', 'bold');
% ylabel('log Signal Amplitude','fontsize',18);
h_legend = legend('SNR = A/\sigma = 0', 'SNR = A/\sigma = 2', 'SNR = A/\sigma = 5', 'Location', 'northeast');
set(h_legend, 'FontSize', 18);