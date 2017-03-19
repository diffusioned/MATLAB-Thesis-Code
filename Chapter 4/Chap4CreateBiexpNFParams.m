% Chap4CreateBiexpNFParams.m
% MATLAB file for plotting data in PhD thesis by Ned Charles available here http://hdl.handle.net/2123/16060
function Chap4CreateBiexpNFParams()

% Create a 70x70 evenly spaced region over the same parameter value range as the original analysis.
NumberOfValues = 4900;
SF1Values = 0:1/69:1;
D1D2RatioValues = 2:18/69:20;
BiexpNFTestParameterArray = zeros(NumberOfValues, 4); % Array parameters are A1, A2, D1, D2
Amp2Values = 1 - SF1Values;
D2Values = 1./D1D2RatioValues;

for i = 1:length(SF1Values)
    for j = 1:length(D1D2RatioValues)
        CurIdx = sub2ind([length(SF1Values) length(D1D2RatioValues)],i,j);
        BiexpNFTestParameterArray(CurIdx,1) = SF1Values(i);
        BiexpNFTestParameterArray(CurIdx,2) = Amp2Values(i);
        BiexpNFTestParameterArray(CurIdx,3) = 1;
        BiexpNFTestParameterArray(CurIdx,4) = D2Values(j);
    end
end

% Save the data
save('YourPath\BiexpNoiseFreeParameters.mat','BiexpNFTestParameterArray');
