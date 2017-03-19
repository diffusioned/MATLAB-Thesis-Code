% Chap2ConvertJacobianIntoMetrics.m
% MATLAB file for fitting simulated data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap2ConvertJacobianIntoMetrics()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Jacobian and create Jacobian condition array as well as JTJ Inv Array
load('YourPath\BiexpFitJacobianArray.mat');
[NFDim SigDim BDim ParamDim] = size(BiexpFitJacobianArray);
JTJInvArray = zeros(NFDim, SigDim,ParamDim,ParamDim);
JacCondArray = zeros(NFDim, SigDim);

for nf = 1:NFDim
    for s = 1:SigDim
        JacCondArray(nf,s) = cond(squeeze(BiexpFitJacobianArray(nf,s,:,:)));
    end
end
save('YourPath\JacCondArray.mat','JacCondArray','-v7.3');

% Also convert to JTJ inverse array
for nf = 1:NFDim
    for s = 1:SigDim
        % Use the QR decomp instead
        RetJac = squeeze(BiexpFitJacobianArray(nf,s,:,:));
        [Q,R] = qr(RetJac,0);
        Rinv = R \ eye(size(R));
        JTJInvArray(nf,s,:,:) = Rinv*Rinv';
    end
end
save('YourPath\JTJInvArray.mat','JTJInvArray','-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Residual Array and convert to RSS
load('YourPath\BiexpResidualArray.mat');
[NFDim SigDim BDim] = size(BiexpResidualArray);
BiexpRSSArray = zeros(NFDim, SigDim);
BiexpRSSArray = sum(BiexpResidualArray.^2,3);
save('YourPath\BiexpRSSArray.mat','BiexpRSSArray','-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the noise sigma associated with each regression
load('YourPath\BiexpRSSArray.mat');
% hard code these, as they aren't in this array
ParamDim = 4;
BDim = 11;
DOFDim = BDim - ParamDim;
RegressionSigmaArray = sqrt(BiexpRSSArray/DOFDim);
save('YourPath\RegressionSigmaArray.mat','RegressionSigmaArray','-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the covariance matrix
load('YourPath\JTJInvArray.mat');
load('YourPath\RegressionSigmaArray.mat');
[NFDim SigDim ParamDim ParamDim] = size(JTJInvArray);
CovMxArray = zeros(NFDim, SigDim, ParamDim, ParamDim);
for nf = 1:NFDim
    for s = 1:SigDim
        CovMxArray(nf,s,:,:) = (RegressionSigmaArray(nf,s).^2).*JTJInvArray(nf,s,:,:); % Need to square the SER array to make variance again
    end
end
save('YourPath\CovMxArray.mat','CovMxArray','-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the correlation matrix
load('YourPath\CovMxArray.mat');
[NFDim SigDim ParamDim ParamDim] = size(CovMxArray);
CorrMxArray = zeros(NFDim, SigDim, ParamDim, ParamDim);
for nf = 1:NFDim
    for s = 1:SigDim
        CurCovMx = squeeze(CovMxArray(nf,s,:,:));
        diagcov = 1 ./ sqrt(diag(CurCovMx));
        DiagCorrMx = spdiags(diagcov,0,4,4);
        CorrMatrix = DiagCorrMx * CurCovMx * DiagCorrMx;
        CorrMxArray(nf,s,:,:) = CorrMatrix(:,:);
    end
end
save('YourPath\CorrMxArray.mat','CorrMxArray','-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('YourPath\CorrMxArray.mat');
[NFDim SigDim ParamDim ParamDim] = size(CorrMxArray);
VIFArray = zeros(NFDim,SigDim,ParamDim);

% Calculate the Variance Inflation Factor(VIF)
for nf = 1:NFDim
    for s = 1:SigDim
        CurCorrMx = squeeze(CorrMxArray(nf,s,:,:));
        VIFArray(nf,s,:) = diag(inv(CurCorrMx));
    end
end
save('YourPath\VIFArray.mat','VIFArray','-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the JTJ condition matrix
load('YourPath\JTJInvArray.mat');
[NFDim SigDim ParamDim ParamDim] = size(JTJInvArray);
JTJConditionArray = zeros(NFDim, SigDim);
for nf = 1:NFDim
    for s = 1:SigDim
        JTJConditionArray(nf,s) = cond(squeeze(JTJInvArray(nf,s,:,:)));
    end
    nf
end
save('YourPath\JTJConditionArray.mat','JTJConditionArray','-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stophere = 1;
