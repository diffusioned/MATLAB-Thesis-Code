# MATLAB-Thesis-Code

This is the code I used to make nearly all of the figures found in my PhD thesis titled "Reliability and Uncertainty in Diffusion MRI Modelling".  This thesis can be downloaded through the University of Sydney's online repository [here](http://hdl.handle.net/2123/16060).  I categorized the figure creation code by the chapters in my thesis and will use this document to link the figures there with the corresponding code.  For many figures in my thesis, many of the figures use the same code so there will be several code files that will cover multiple figures.

### Chapter 1 - Introduction ###

**Figure 1** was created using the drawing tools in Microsoft Office suite.

**Figure 2** - I wasn't able to locate the code for this figure and it was probably buried deep in the folder structure of my school laptop.  I'm sure you are well aware how your file and folder naming conventions change over the course of your work, academia or industrial.  You can recreate this file by following the explanation in the figure caption to add noise to a line and fit a linear model to it.

[Chap1ShowRicianPDFs.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%201/Chap1ShowRicianPDFs.m) is the file to create **Figure 3**.

**Figure 4** was created with a little bit of Photoshop magic.  First, I plotted the dashed line with a log plot.  Then I did another scatter plot of about 10,000 red points with added magnitude noise, demonstrating the Rician signal bias.  I then copied the scatter plot and overlaid it on top of the dashed line plot and adjusted the transparency.  Adding magnitude noise to a signal will be demonstrated in Chapter 2.

[Chap1OverFitCurveOnLine.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%201/Chap1OverFitCurveOnLine.m) is the file to create **Figure 5**.


### Chapter 2 - The Biexponential Model ###

**Figure 6** was a straight forward scatter plot based on the values in Table 1, color coded by tissue type.

**Figure 7 and 8** were both created in the same process as Figure 4.

[Chap2LinePlotOfSignalRange.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2LinePlotOfSignalRange.m) was used to create the line plot in **Figure 9**.  The blue shaded area was added via Photoshop like Figure 4.

For the rest of this chapter, as well as Chapters 3 and 4, simulated noise-free data will be needed.  For Chapter 2, however, this data consisted of ten million individual signals from simulated noisy measurements.  This produces an 800 MB .mat file, which I won't upload to GitHub.  There will be even larger files being created later on when saving the diagnostic info from the regression fits.  You can probably much fewer signals if you want test this out.  In Chapter 4, I reduced the number of noise-free signals to 4900 and you can compare the results there to Chapter 2, which are very similar.

[Chap2CreateNoiseFreeParameters.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2CreateNoiseFreeParameters.m) is the first file I used to create my test set.  In it I created 50,000 random unique parameter combinations for a biexponential signal.  Set the number of noise-free signals you'd like to create here.  

[Chap2CreateNoiseFreeSignals.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2CreateNoiseFreeSignals.m) then takes each of these random parameter combinations and creates a signal for it based on the simulated diffusion weights or *b*-values you provide within.  Adding more diffusion weightings will create more data, but you can adjust the values to whatever you want.

[Chap2CreateNoisySignals.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2CreateNoisySignals.m) then takes each of these noise-free signals and creates multiple noisy signals based on the number you set within.  The simulated added noise is based on a magnitude signal, giving the data a Rician signal bias.
Now that you have these data you can do a NLLS regression fit of them to compare the fitted parameter values to the original noise-free parameter values.

[Chap2NLLSBiexpFit.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2NLLSBiexpFit.m) is the file that I used to fit both biexponential and monoexponential models to my simulated noisy data.  If you are fitting ten million noisy signals, be forewarned that this file will take a LONG time to complete.  I have a parfor statement in this file so that MATLAB will speed up the analysis by running the regression fitting loop in parallel.  However, when running this file with 12 parallel processors performed on the University of Sydney's High Performance Cluster, it still took about 30 hours.  Note, this amount of signals will also produce a saved Jacobian array of about 3-4 GB.  This size of file will require 64-bit MATLAB to be installed.  After running this fitting file, the rest of the analysis is much faster.

[Chap2ConvertJacobianIntoMetrics.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2ConvertJacobianIntoMetrics.m) As I explain in my [Nonlinear Regression Primer](https://github.com/diffusioned/diffusioned.github.io/blob/master/_posts/2016-09-05-Nonlinear_Regression.md), the Jacobian matrix, returned in the above biexponential fitting file, is very useful in producing diagnostic information about the fitting parameters.  As I stated above, the array of Jacobian matrices is very large, to this file uses these matrices to create and save various diagnostic measures.  This will also create a lot of data, and you can comment out portions of that you are interested in.

[Chap2DisplayBiexpFitParameterMaps.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2DisplayBiexpFitParameterMaps.m) is the file that displayed the parameter fit maps for the biexponential and monoexponential fit parameters estimated in the fitting function.  This file was used to create the nine individual parameter pseudocolor plots in **Figure 10** of my thesis, which were then stitched together in Photoshop.  A similar process was used for the robust measures plotted with this file to create **Figure 13**.  This file also produced plot of the monoexponential parameter estimates in the top row of **Figure 14** as well as **Figure 15**.

[Chap2PlotHistograms.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2PlotHistograms.m) is the file used to plot the histograms in **Figures 11 and 12**.  This file also contains the optional sections to plot the same histograms with confidence intervals added in **Figures 30 and 31**

[Chap2PlotFiveSignals.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2PlotFiveSignals.m) creates the line plot at the bottom of **Figure 14**.

[Chap2DisplayMonoexpSER.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2DisplayMonoexpSER.m) was used to create **Figure 16**.

[Chap2DisplayBiexpFitDiagnostics.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2DisplayBiexpFitDiagnostics.m) displays a lot of the diagnostic measures that were created earlier through the Jacobian matrix calculations.  This covers the biexponential SER plot and histogram in **Figure 17**, the QQ plots in **Figure 18**, the Jacobian condition plots in **Figures 24 and 25**, the mean and median parameter standard deviation plots in **Figures 26 and 27**, and the confidence interval percentage plots in **Figure 28**.

[Chap2NLLSBiexp2SDFit.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2NLLSBiexp2SDFit.m) is a modified version of *Chap2NLLSBiexpFit.m* above, which incorporated code to fit only signal points that are greater than two times the standard deviation of the simulated noise.

[Chap2DisplaySignalAverageSNR.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2DisplaySignalAverageSNR.m) created the signal-averaged SNR plots in **Figures 19 and 20**.

[Chap2Display2SDMaxSignalValue.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2Display2SDMaxSignalValue.m) took the data from the 2SD fit file above and plotted the maximum signal value for all signals tested.

[Chap2CompareFulland2SDFitOnMonoSER.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2CompareFulland2SDFitOnMonoSER.m) also used the 2SD fit data to create **Figure 22**.

**Figure 23** was created using the same code as in *Chap2DisplayBiexpFitParameterMaps.m* except with the 2SD fit data

[Chap2Plot95PctGaussianConfidenceInterval.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2Plot95PctGaussianConfidenceInterval.m) plotted the Gaussian distribution seen in **Figure 29**.

The bootstrap test set in Chapter 2 consisted of only 400 noise-free signals.  To replicate these, modify the *Chap2CreateNoiseFreeSignals.m* file above with these lines instead:

	NumberOfValues = 400;
	SF1Values = 0.025:0.05:0.975;
	D1D2RatioValues = 2.45:0.9:19.55;

Then create the noise-free signals using the *Chap2CreateNoiseFreeSignals.m* file and then the simulated noisy signal measurements using *Chap2CreateNoisySignals.m*.  Note I only created 25 noisy measurements for each noise-free signal in the bootstrap section.  

[Chap2NLLSBootstrapBiexpFit.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2NLLSBootstrapBiexpFit.m) is the file that was used to fit these noisy signal measurements.  This file created and fit 1000 bootstrap samples of each  noisy measurement regression fit.

[Chap2DisplayBootstrapBiexpFitMaps.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2DisplayBootstrapBiexpFitMaps.m) is the file that created **Figures 32, 33, 34, and 35**.

[Chap2CompareRSSFitContours.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%202/Chap2CompareRSSFitContours.m) displayed the RSS contour values seen in **Figures 36 and 37**.  You can uncomment the various sections to produce the plot you want or tinker with the values to display your own.

Finally, **Figure 38** was produced using the contours from one plot in Figure 10 overlaid on Figure 6.

### Chapter 3 - The Kurtosis Model ###

[Chap3PlotMonoexpSignalWithKurtosisParameters.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%203/Chap3PlotMonoexpSignalWithKurtosisParameters.m) plotted **Figure 39**.

The data set for the biexponential signals fit with the kurtosis model was the same as used in Chapter 2.  This chapter had an additional monoexponential test set.

[Chap3CreateMonoexpNoiseFreeParameters.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%203/Chap3CreateMonoexpNoiseFreeParameters.m) was the function that created the monoexponential test parameters.  Use the noise-free signal and noisy signal .m files from Chapter 2 to create the measurement data.  I created 1000 noise-free signals and then tested 200 noisy signal measurements for each of those signals.

[Chap3NLLSKurtosisFitToBiexpData.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%203/Chap3NLLSKurtosisFitToBiexpData.m) and [Chap3NLLSKurtosisFitToMonoexpData.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%203/Chap3NLLSKurtosisFitToMonoexpData.m) were the two files that fit the kurtosis model to the two data sets.

[Chap3DisplayKurtosisFitParamsOnBiexpTruth.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%203/Chap3DisplayKurtosisFitParamsOnBiexpTruth.m) was the file that displayed the fitted parameter values to the biexponential test set to create **Figure 40**.

[Chap3DisplayKurtosisFitJacobian.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%203/Chap3DisplayKurtosisFitJacobian.m) was the file that displayed the Jacobian condition plot in **Figure 41**.  Note that the condition number has to be calculated via the same process as in *Chap2ConvertJacobianIntoMetrics.m*.

**Figure 42** was displayed using *Chap3DisplayKurtosisFitParamsOnBiexpTruth.m* for the higher SNR fit data.

[Chap3KurtHistLinePlotThreeSignals.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%203/Chap3KurtHistLinePlotThreeSignals.m) created the top left line plot in **Figure 43**.  The top right plot is simply Figure 41.

[Chap3ParameterHistogramDisplay.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%203/Chap3ParameterHistogramDisplay.m) created the histograms in Figure 43.

The bootstrap signals were the same limited test set from the original biexponential test set using the code snippet above for Chapter 2.

[Chap3DisplayBootstrapHistograms.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%203/Chap3DisplayBootstrapHistograms.m) displayed the bootstrap histograms shown in **Figure 44**.

The > 2SD fitting was performed in Chapter 3 using the same modifications to the fit code as *Chap2NLLSBiexp2SDFit.m*

[Chap3CompareFulland2SDFitOnKurtSER.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%203/Chap3CompareFulland2SDFitOnKurtSER.m) used the 2SD data to produce **Figure 45**.

[Chap3DisplayKurtFitOnBiexpSER.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%203/Chap3DisplayKurtFitOnBiexpSER.m) produced **Figure 46**.

[Chap3DisplayKurtosisFitParamsOnMonoexpTruth.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%203/Chap3DisplayKurtosisFitParamsOnMonoexpTruth.m) produced the final figure in the chapter, **Figure 47**.

### Chapter 4 - Model Selection ###

The model selection has a reduced test set for the noise-free biexponential signals.

[Chap4CreateBiexpNFParams.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%204/Chap4CreateBiexpNFParams.m) produces the noise-free parameter combinations, and then the noise-free signal and noisy signal measurement creation code from Chapter 2 can be used again.  However, Chapter 4 has two sets of diffusion weightings being tested.  Add this code to *Chap2CreateNoiseFreeSignals.m* to set the values. 

	% BValueArray = [0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4, 12.8, 25.6]; % Chapter 4 11 b-values
	BValueArray = [0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6]; % Chapter 4 7 b-values

The monoexponential signals are the same as Chapter 3 with the modification of this same code to create two data sets, as well.  The fitting code is similar to the previous chapters with the addition of code that calculates the leave-one-out cross-validation (LOOCV).

[Chap4NLLSandLOOCVFit3ModelsToBiexpTruth.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%204/Chap4NLLSandLOOCVFit3ModelsToBiexpTruth.m) was a sample of that code used for the different biexponential data sets.  It can also be used on the monoexponential data.

The first several plots in the model selection chapter use the same code for the different data sets.  This code produces nine individual plots and then Photoshop was used to stitch the plots together into one image.

[Chap4DisplaySelectionPct11BValuesSNR25](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%204/Chap4DisplaySelectionPct11BValuesSNR25Biexp.m) was the code that produces the nine individual plots shown for 11 diffusion weightings and SNR 25 in **Figure 48**.  Similar plots are created for the fits of 11 diffusion weightings and SNR 200 for **Figure 49**, 7 diffusion weightings and SNR 25 for **Figure 50**, and 7 diffusion weightings and SNR 200 for **Figure 51**.

[Chap4DisplaySelectionPct11BValuesSNR25Monoexp.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%204/Chap4DisplaySelectionPct11BValuesSNR25Monoexp.m) was a similar example for the monoexponential histogram plots.  This file produced **Figure 52** for 11 diffusion weightings and SNR 25, with similar plots for 11 diffusion weightings and SNR 200 in **Figure 53**,  7 diffusion weightings and SNR 25 in **Figure 54**, and 7 diffusion weightings and SNR 200 in **Figure 55**.

[Chap4DisplayDeltaAIC11BVSNR25.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%204/Chap4DisplayDeltaAIC11BVSNR25.m) was the file that produces the ΔAIC plots in **Figure 56**.  This section can also be used to create the plot for the SNR 200 data in **Figure 59**.  The second portion of this file creates the plots displaying the difference between ΔAIC and ΔAICc in **Figure 62**.

[Chap4DisplayDeltaAICCaseStudies11BVSNR25.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%204/Chap4DisplayDeltaAICCaseStudies11BVSNR25.m) produced the plots in **Figures 57 and 58**.  It can also be used to create the plots for the SNR 200 data in **Figures 60 and 61**

[Chap4DisplaySelectionPctAICVsFTest.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%204/Chap4DisplaySelectionPctAICVsFTest.m) created **Figure 63**.

[Chap4DisplaySelectionPctVsSkewedParameters.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%204/Chap4DisplaySelectionPctVsSkewedParameters.m) was the last file that created **Figures 64, 65, 66, and 67**.

### Chapter 5 - Actual Tissue Testing ###

Unfortunately, the data set I used in my thesis is not mine to share, so you won't be able to reproduce the plots seen there.  However, if you are that interested in this work, you probably have your own data.  There are also freely available data sets out on the internet to work with.  I have produced a few code files here that demonstrate how I displayed my data in my thesis.

All data was fit with the same code as *Chap2NLLSBootstrapBiexpFit.m* for the biexponential model, and this code can be used for the kurtosis and monoexponential fits with some modifications to the array sizes and model equations.

[Chap5DisplayAverageSNRForTissueData.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%205/Chap5DisplayAverageSNRForTissueData.m) was the file that produced **Figure 69**.

[Chap5DisplayModelFitToTissueData.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%205/Chap5DisplayHistogramsForFitTissueData.m) was the file that produced **Figures 70 and 71** in my thesis.  

[Chap5DisplayHistogramsForFitTissueData.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%205/Chap5DisplayModelFitToTissueData.m) is an example of how the parameter estimates were displayed as a histogram, specifically the biexponential parameter results in **Figure 74**.  It can also be modified to display the monoexponential and kurtosis parameter estimates in **Figures 72 and 73**.

[Chap5DisplayBiexpFitROIBootCI.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%205/Chap5DisplayBiexpFitROIBootCI.m) was an example of how to create the line plots displaying the initial estimates as a dot with the bootstrap range overlapped as a line.  This particular file was used to create **Figure 80**, but it can also be modified to produce **Figures 75, 76, 77, 78, and 79**.

**Figure 81** can be replicated by plotting the original noisy measurement along with a selected bootstrap sample measurement.

[Chap5DisplayBiexpFitWithNormalityBootElim.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%205/Chap5DisplayBiexpFitWithNormalityBootElim.m) produced the biexponential results in  **Figures 82 and 83**, and can be modified to produce the kurtosis results in **Figure 84**.

**Figure 85** can be replcated using *Chap5DisplayAverageSNRForTissueData.m* with the selected best model from the your AIC results.

[Chap5DisplayBiexpFitWithEliminationROIDelAIC10.m](https://github.com/diffusioned/MATLAB-Thesis-Code/blob/master/Chapter%205/Chap5DisplayBiexpFitWithEliminationROIDelAIC10.m) produced the results where biexponential fits with ΔAIC < 10 were removed, namely **Figures 86, 88, and 89**.  It can be modified to produced the kurtosis fits results in **Figures 87, 90, and 91**.

### I hope some of this code helps you in your work.###