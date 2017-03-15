# MATLAB-Thesis-Code

This is the code I used to make nearly all of the figures found in my PhD thesis titled "Reliability and Uncertainty in Diffusion MRI Modelling".  This thesis can be downloaded through the University of Sydney's online repository [here](http://hdl.handle.net/2123/16060).  I categorized the figure creation code by the chapters in my thesis and will use this document to link the figures there with the corresponding code.  For many figures in my thesis, many of the figures use the same code so there will be several code files that will cover multiple figures.

### Chapter 1 - Introduction ###

**Figure 1** was created using the drawing tools in Microsoft Office suite.

**Figure 2** - I wasn't able to locate the code for this figure and it was probably buried deep in the folder structure of my school laptop.  I'm sure you are well aware how your file and folder naming conventions change over the course of your work, academia or industrial.  You can recreate this file by followin the explanation in the figure caption to add noise to a line and fit a linear model to it.

[Chap1ShowRicianPDFs.m](link) is the file to create **Figure 3**.

**Figure 4** was created with a little bit of Photoshop magic.  First, I plotted the dashed line with a log plot.  Then I did another scatter plot of about 10,000  red points with added magnitude noise, demonstrating the Rician signal bias.  I then copied the scatter plot and overlaid it on top of the dashed line plot and adjusted the transparency.  Adding magnitude noise to a signal will be demonstrated in Chapter 2.

[Chap1OverFitCurveOnLine.m](link) is the file to create **Figure 5**.


### Chapter 2 - The Biexponential Model ###

**Figure 6** is a straight forward scatter plot based on the values in Table 1, color coded by tissue type.

**Figure 7 and 8** were both created in the same process as Figure 4.

[Chap2LinePlotOfSignalRange.m](link) was used to create the line plot in **Figure 9**.  The blue shaded area was added via Photoshop like Figure 4.

For the rest of this chapter, as well as Chapters 3 and 4, simulated noise-free data will be needed.  For Chapter 2, however, this data consisted of ten million individual signals from simulated noisy measurements.  This produces an 800 MB .mat file, which I won't upload to GitHub.  There will be even larger files being created later on when saving the diagnostic info from the regression fits.  You can probably much fewer signals if you want test this out.  In Chapter 4, I reduced the number of noise-free signals to 4900 and you can compare the results there to Chapter 2, which are very similar.

[Chap2CreateNoiseFreeParameters.m](link) is the first file I used to create my test set.  In it I create 50,000 random unique parameter combinations for a biexponential signal.  Set the number of noise-free signals you'd like to create here.  

[Chap2CreateNoiseFreeSignals.m](link) then takes each of these random parameter combinations and creates a signal for it based on the simulated diffusion weights or *b*-values you provide within.  Adding more diffusion weightings will create more data, but you can adjust the values to whatever you want.

[Chap2CreateNoisySignals.m](link) then takes each of these noise-free signals and creates multiple noisy signals based on the number you set within.  The simulated added noise is based on a magnitude signal, giving the data a Rician signal bias.

Now that you have these data you can do a NLLS regression fit of them to compare the fitted parameter values to the original noise-free parameter values.