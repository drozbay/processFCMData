%% 1D DeltaF/F %%
% Performs the dF/F calculation on a fluorescence transient signal using
% the method described in Jia et al. 2011. Nature Protocols
% Inputs:
% tIn - Vector of evenly-spaced time points corresponding to fluorescence
% measurements
% fIn - Vector of fluorescence measurements for each time point
% t0 - Time constant of exponentially-weighted noise filter (typical: 0.2)
% t1 - Width of minimum sampling time window (typical: 0.75)
% t2 - Temporal lag of time window (typical: 2)
% Outputs:
% dFOut - Resulting deltaF/F signal (same size as fIn)
% f0 - Windowed minimum fluorescence values (same size as fIn)
% f1 - Windowed mean standard deviation (same size as fIn)

function [dFOut,f0,f1] = deltaFCalc(tIn, fIn,t0,t1,t2)
if size(fIn,1)==1
    fIn=fIn';
    rememberToTranspose = 1;
else
    rememberToTranspose = 0;
end
dT = tIn(2)-tIn(1);
n0 = t0/dT;
n1 = round(t1/dT);
n2 = round(t2/dT);

fInPad = padarray(fIn,n2+n1,'replicate','both');
f0 = zeros(length(fIn),1);
f1 = zeros(length(fIn),1);
meanF = zeros(n2,1);
stdF = zeros(n2,1);

for ii = 1:length(tIn)
    xx = ii+n1:ii+n1+n2;
    for nn = 1:n2
        meanF(nn) = mean(fInPad(xx(nn)-n1:xx(nn)+n1));
        stdF(nn) = std(fInPad(xx(nn)-n1:xx(nn)+n1),1);
    end
    f0(ii) = min(meanF);
    f1(ii) = mean(stdF);
end
f1(f1<mean(f1)/100) = mean(f1);
fR = (fIn - f0)./f0;
% fR = ((fIn - f0)./f0)/mean(f1);

alpha = exp(-n0);
dFOut = filter(alpha,[1 alpha-1], fR);
if rememberToTranspose
    dFOut = dFOut';
end