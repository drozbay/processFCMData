close all
clear all

minSize=20;  % Minimum area threshold
dFlim=2; % Maximum delta F (used to exclude artifacts)
tWin=5; % Time window
sdThreshold=5; % Criterion for response
pWin=20; % Number of sample to include before and after peak for windowed plot


%dt_sample=0.545; %Seconds between samples, good for 5
%dt_sample=0.404; %Seconds between samples, good for 11
%dt_sample=0.520; %Seconds between samples, good for 13
dt_sample=0.303; %Seconds between samples, good for 15
% dt_sample=0.848; %Seconds between samples, good for 17

% Get image path and filename from user
[fname,pname,nCancel] = uigetfile({'*.tif;*.tiff'},'Select the TIMELAPSE file...');
if nCancel
    inputPath = [pname,fname];
    pnameStart = pname;
else
    error('Cancelled')
end

%% Import images
% Get image information
infoInput = imfinfo(inputPath);
% Number of frames to drop
frameDrop = 200;
% numImages is final frame to collect
numImages = length(infoInput)-frameDrop;
xInput = infoInput.Width;
yInput = infoInput.Height;
% Set up matrix
imInputDb = double(zeros(yInput,xInput,numImages));
xInputSc = size(imInputDb,2);
yInputSc = size(imInputDb,1);
% Read images into stack
hWait = waitbar(0, sprintf('Importing %d images...',numImages));
for ii = 1:numImages
    imInputDb(:,:,ii) = imread(inputPath,'tif',ii);
    waitbar(ii/numImages,hWait);
end
close(hWait);

%% Get responsiveness ROI image path and filename from user
BW=imread([inputPath(1:end-4) 'ROI.tif']);
% Create distance to center map
D = -bwdist(~BW);
% Watershed
Ld = watershed(D);
% Use the watershed lines to create borders
BW2 = BW;
BW2(Ld==0) = 0;
% Remove all rois smaller than area threshold
% Fill holes inside ROIS
BW3 = imfill(bwareaopen(BW2, minSize),'holes');
cc = bwconncomp(BW3);
% Make color for plot
coloredLabels = label2rgb(bwlabel(BW3), 'hsv', 'k', 'shuffle');

% Plot region masks
figure(1);
subplot(1,3,1);
imshow(BW);
title('Response threshold >a*SD');
% Watershed
subplot(1,3,2);
imshow(BW3);
title('Watershed ROI');
% Color plot
subplot(1,3,3);
imshow(coloredLabels);
title('Individual areas');

%% Extract all areas from above
% Data structure for ROIs
fullROIData = struct(...
    'ROINum',[],...
    'intensity',zeros(1,numImages),...
    'deltaF',zeros(1,numImages),...
    'stdVal',[],...
    'medianVal',[],...
    'hasResponse',zeros(1,numImages,'logical'),...
    'responseIdx',zeros(1,numImages),...
    'deltaFPre',[],...
    'deltaFPeak',[]);

% Get get mean intensity of each ROI for entire sequence
tempCurrentImage=zeros(xInputSc,yInputSc);
numROIs=length(regionprops(cc, imInputDb(:,:,1), 'all'));
hWait = waitbar(0, sprintf('Processing %d images...',numImages));
for ii = 1:numImages
    tempCurrentImage = imInputDb(:,:,ii);
    ROIprop = regionprops(BW3, tempCurrentImage, 'all');
    for jj=1:numROIs
        fullROIData(jj).intensity(ii)=ROIprop(jj).MeanIntensity;
    end
    waitbar(ii/numImages,hWait);
end
% Assign a number for each ROI
for jj = 1:numROIs
    fullROIData(jj).ROINum = jj;
end
clear temp*;
close(hWait);

%% Processing of data
timeVector=linspace(0,(numImages-1)*dt_sample,numImages);
% DeltaF/F Calculation
t0 = 0.5;
t1 = dt_sample*3;
t2 = dt_sample*6;
hWait = waitbar(0, sprintf('Performing deltaF/F on %d ROIs...',numROIs));
for jj = 1:numROIs
    fullROIData(jj).deltaF = deltaFCalc(timeVector,...
        fullROIData(jj).intensity,t0,t1,t2);
    waitbar(jj/numROIs,hWait);
end
close(hWait);
%%
jj=1;
[fullROIData(jj).deltaF,f0,f1] = deltaFCalc(timeVector,...
    fullROIData(jj).intensity,t0,t1,t2);
figure(2); clf;
subplot(2,1,1); hold on;
plot(timeVector,fullROIData(jj).intensity,'-k') 
plot(timeVector,f0,'-r');
xlim([0 timeVector(end)]);
ylim([min(fullROIData(jj).intensity),...
    max(fullROIData(jj).intensity)]);
xlabel('seconds');

subplot(2,1,2);
plot(timeVector,fullROIData(jj).deltaF,'-b');
xlim([0 timeVector(end)]);
ylim([min(fullROIData(jj).deltaF) 0.5]);
ylabel('deltaF/F');
xlabel('seconds');

%% Get an estimate of the median standard deviation
for jj=1:numROIs
    fullROIData(jj).stdVal = std(fullROIData(jj).deltaF,0);
    fullROIData(jj).medianVal = median(fullROIData(jj).deltaF);
end
clear tempStdVals;
figure(3); clf;
histogram([fullROIData.stdVal]);
xlabel('Standard deviation');
title('Distribution of standard deviation for traces');

%% Find responses using standard deviation threshold
isResponse=zeros(numImages,numROIs);
for jj=1:numROIs
    tempStdVal = fullROIData(jj).stdVal;
    tempMedVal = fullROIData(jj).medianVal;
    %Calculate responses
    for ii=1:numImages
        if fullROIData(jj).deltaF(ii) > ...
                sdThreshold*tempStdVal+tempMedVal
            isResponse(ii,jj)=1;
        end
    end
    % Remove repeated vaues from isResponse vector
    tempCleanResponses = [0; diff(isResponse(:,jj))];
    tempCleanResponses(tempCleanResponses<0) = 0;
    isResponse(:,jj) = isResponse(:,jj).*tempCleanResponses;
    % Remove responses that occur within tWin from a previous one
    for ii = 1:numImages-tWin
       if isResponse(ii,jj)
           isResponse(ii+1:ii+tWin,jj) = 0;
       end
    end
end
clear temp*;

fprintf('Total responses: %d\n',numel(isResponse(isResponse==1)));

% Display list of ROIs and indicate responses
figure(4); 
set(gcf,'Units','Inches','Position',[5 0.5 5 9]);
clf;
hold on;
traceSep = 0.5;
for jj=1:numROIs
    tempStdVal = fullROIData(jj).stdVal;
    tempMedVal = fullROIData(jj).medianVal;
    plot(timeVector, fullROIData(jj).deltaF - tempMedVal+(jj-1)*traceSep);
    tempIsResponse = isResponse(:,jj);
    tempIsResponse(tempIsResponse==0)=10e6;
    plot(timeVector,tempIsResponse/4+(jj-1)*traceSep,'.r','markersize',10);
    plot(timeVector,repmat(tempStdVal,numImages,1)+(jj-1)*traceSep,...
        '--b','linewidth',0.1);
    plot(timeVector,zeros(1,numImages)+(jj-1)*traceSep,...
        '-b','linewidth',0.1);
    text(-timeVector(end)/100,(jj-1)*traceSep,sprintf('%d',jj),...
        'HorizontalAlignment','Right');
end
xlim([0 timeVector(end)]);
ylim([-traceSep (numROIs+1)*traceSep]);
set(gca,'ytick',[]);
clear temp*;

%% Count the number of responses
for jj=1:numROIs
    if sum(isResponse(:,jj))>=1
        fullROIData(jj).hasResponse = 1;
        fullROIData(jj).responseIdx = find(isResponse(:,jj));
    else
        fullROIData(jj).hasResponse = 0;
    end
end
% Remove ROIs that do not have a response
tempIdx = find([fullROIData.hasResponse]==1);
respROIData = fullROIData(tempIdx);
clear temp*;

%% Get pre-transient and peak-transient values
deltaFTraces = zeros(length([fullROIData.hasResponse]==1),pWin*2+1);
tempTraceIdx = 0;
for jj=1:length(respROIData)
    tempNumResponses = length(respROIData(jj).responseIdx);
    tempDeltaFPre = zeros(tempNumResponses,1);
    tempDeltaFPeaks = zeros(tempNumResponses,1);
    tempDeltaF = padarray(respROIData(jj).deltaF,[tWin 0],'replicate');
    for ii=1:tempNumResponses
        tempIdx = respROIData(jj).responseIdx(ii)+tWin;
        % Get the peak value tWin around the response
        tempDeltaFPeaks(ii) = ...
            max(tempDeltaF(tempIdx-tWin:tempIdx+tWin));
        % Get the mean local value tWin ahead of the response
        tempDeltaFPre(ii) = ...
            mean(tempDeltaF(tempIdx-tWin-1:tempIdx-tWin+1));
        % Get windowed trace of transient
        if (tempIdx-tWin-pWin>0) && (tempIdx-tWin+pWin<numImages)
            tempTraceIdx = tempTraceIdx + 1;
            deltaFTraces(tempTraceIdx,1:pWin*2+1) = ...
                tempDeltaF(tempIdx-pWin:tempIdx+pWin);
        end
    end
    deltaFTraces = deltaFTraces(1:tempTraceIdx,:);
    respROIData(jj).deltaFPeak = tempDeltaFPeaks;
    respROIData(jj).deltaFPre = tempDeltaFPre;
end
clear temp*;
% Plot histogram of all transients
figure(5); clf; hold on;
h1 = histogram(vertcat(respROIData.deltaFPeak));
h2 = histogram(vertcat(respROIData.deltaFPre));
h1.BinWidth = h2.BinWidth;
title('Summary of all transients');
xlabel('deltaF/F');
legend('Peak Transient Values','Pre-Transient Values');

%% Plot summary
roiCMap = hsv(length(respROIData))*0.8;
figure(6); clf;
subplot(2,2,1); cla;
BWwResponse = ismember(labelmatrix(cc), [respROIData.ROINum]);
BWwResponseRGB = label2rgb(bwlabel(BWwResponse), roiCMap, 'k');
imshow(BWwResponseRGB,'initialmagnification','fit'); hold on;
for jj = 1:length(respROIData)
    tempIdx =respROIData(jj).ROINum;
    text(ROIprop(tempIdx).Centroid(1),ROIprop(tempIdx).Centroid(2),...
        sprintf('%d',tempIdx),'Color','w');
end

subplot(2,2,2); cla; hold on;
traceSep = 0.4;
for jj=1:length(respROIData)
    tempStdVal = respROIData(jj).stdVal;
    tempMedVal = respROIData(jj).medianVal;
    plot(timeVector,respROIData(jj).deltaF-tempMedVal+(jj-1)*traceSep,...
        '-','Color',roiCMap(jj,:))
    text(-timeVector(end)/100,(jj-1)*traceSep,...
        sprintf('%d',respROIData(jj).ROINum),...
        'Color',roiCMap(jj,:),'HorizontalAlignment','Right');
    % Add arrows at responses
    for ii=1:length(respROIData(jj).responseIdx)
        tempTimePoint = respROIData(jj).responseIdx(ii)*dt_sample;
        tempH = annotation('arrow');
        set(tempH,'parent', gca, ...
            'position', [tempTimePoint (jj-1.5)*traceSep 0 traceSep*0.5], ...
            'HeadLength', 2, 'HeadWidth', 2, 'HeadStyle', 'plain', ...
            'Color', roiCMap(jj,:)*0.6, 'LineWidth', 0.5);

    end
end

xlim([0 timeVector(end)]);
ylim([-traceSep (length(respROIData)+1)*traceSep]);
set(gca,'YTick',[]);
xlabel('Time (sec)');
ylabel(['Region #',newline]);
grid on; grid minor;

% Show the spike timecourses
subplot(2,2,3);
hold on
plot((1:pWin*2+1)*dt_sample,deltaFTraces,'-')
plot([pWin+1 pWin+1]*dt_sample,[min(deltaFTraces(:)) max(deltaFTraces(:))],'-r')
xlabel('Time (sec)');
title('Ca2+ transient timecourses');
xlim([1 pWin*2+1]*dt_sample);
xlabel('seconds');
ylabel('DeltaF/F');
grid on;

% Show spikes raster plot
subplot(2,2,4); cla;
spikeTimes = vertcat(respROIData.responseIdx)*dt_sample;
tempDeltaFPeaks = vertcat(respROIData.deltaFPeak);
scatter(spikeTimes,tempDeltaFPeaks,20,'s','filled');
xlim([0,timeVector(end)]);
title('Calcium Transients');
xlabel('Time (sec)');
ylabel('DeltaF/F');
grid on; grid minor;

clear temp*;

%% Show the spike timecourses, shifted
figure(7); clf;
maxdFFs=max(deltaFTraces(:));
mindFFs=min(deltaFTraces(:));
hold on;
for ii=1:size(deltaFTraces,1)
    plot((1:pWin+pWin+1)*dt_sample,deltaFTraces(ii,:)+ii*0.5*(maxdFFs-mindFFs),'-')
end
plot([pWin+1 pWin+1]*dt_sample,[min(deltaFTraces(:)) max(deltaFTraces(:))],'-r')
xlim([0 pWin+pWin+1]*dt_sample);

%% Get the max projection images of time points with responses
imInput_Filt = zeros(size(imInputDb));
imInput_Filt_NoResp = imInput_Filt;

% Filter images
for ii=1:numImages
    imInput_Filt(:,:,ii) = imgaussfilt(medfilt2(imInputDb(:,:,ii),[1,1]),1);
end

% Subtract background
% Remove images in window around transients
imInput_Filt_noResp = padarray(imInput_Filt,[0 0 tWin],0);
numTransients = length(vertcat(respROIData.responseIdx));
tempRespIdx_tWin = sort(vertcat(respROIData.responseIdx))+tWin;
for jj=1:numTransients
    imInput_Filt_noResp(:,:,tempRespIdx_tWin(jj)-tWin:tempRespIdx_tWin(jj)+tWin) =...
        zeros(size(imInput_Filt,1),size(imInput_Filt,2),tWin*2+1);
end
%%
imInput_Filt_noResp = imInput_Filt_noResp(tWin+1:end-tWin);
noRespRemoveIdx = squeeze(max(max(imInput_Filt_noResp,[],1),[],2))>0;
imInput_Filt_noResp_Short = imInput_Filt_noResp(:,:,~noRespRemoveIdx);
imInput_Filt_wResp_Short = imInput_Filt_noResp(:,:,noRespRemoveIdx);
%%
tempRespIdx = sort(vertcat(respROIData.responseIdx));

imMeanFull = mean(imInput_Filt,3);
imPeakTransients = zeros(xInput,yInput,numTransients);
figure(9); set(gcf,'Units','Inches','Position',[6 3 5.5 5]); clf;
for ii=1:size(deltaFTraces,1)
    imPeakTransients(:,:,ii) = max(imInput_Filt(:,:,tempRespIdx(ii)-1:tempRespIdx(ii)+1),[],3)-imMeanFull;
end
imMaxProjPeaks = max(imPeakTransients,[],3);
imagesc(imMaxProjPeaks);
colormap(gray);
title('Max projection');
imwrite(uint16(imMaxProjPeaks./max(imMaxProjPeaks).*2^16-1),'imMaxProjPeaks.tif');

imInput_Filt_noResp_Short = uint16(imInput_Filt_noResp_Short./max(imInput_Filt_noResp_Short(:)).*2^16-1);
imwrite(imInput_Filt_noResp_Short(:,:,1),'imInput_Filt_NoResp.tif');
for ii = 2:size(imInput_Filt_noResp_Short,3)
    imwrite(imInput_Filt_noResp_Short(:,:,ii),'imInput_Filt_NoResp.tif','writemode','append');
end

imInput_Filt_wResp_Short = uint16(imInput_Filt_wResp_Short./max(imInput_Filt_wResp_Short(:)).*2^16-1);
imwrite(imInput_Filt_wResp_Short(:,:,1),'imInput_Filt_wResp.tif');
for ii = 2:size(imInput_Filt_wResp_Short,3)
    imwrite(imInput_Filt_wResp_Short(:,:,ii),'imInput_Filt_wResp.tif','writemode','append');
end


%%
%Get the fwhm
fwhm_per_spike=[];
no_s_fwhm=0;
for ii=1:no_spikes
    
    %The fwhm in MATLAB exchange does not work with a few points
    %Becaus of that problem D. Restrepo wrote this code
    %Please not this code is hard coded and will only work for other
    %transients if the number of points in the transient is adjusted
    %accordingly
    
    %Get fwhm
    try
        %Get the baseline
        dFFo=mean(deltaFTraces(ii,ii_pre+1-6:ii_pre+1-2));
        
        
        %Find the maximum
        [maxdFFo maxii]=max(deltaFTraces(ii,ii_pre+1-6:ii_pre+10));
        
        %Find the half point up
        jj=1;
        while deltaFTraces(ii,ii_pre+1-6+jj-1)<dFFo+0.5*(maxdFFo-dFFo)
            jj=jj+1;
        end
        slope1=(deltaFTraces(ii,ii_pre+1-6+jj-1)-deltaFTraces(ii,ii_pre+1-6+jj-2))/dt_sample;
        tbef=(jj-1)*dt_sample+((dFFo+0.5*(maxdFFo-dFFo)-deltaFTraces(ii,ii_pre+1-6+jj-2))/slope1);
        
        %Find half point down
        jj=maxii;
        while deltaFTraces(ii,ii_pre+1-6+jj-1)>dFFo+0.5*(maxdFFo-dFFo)
            jj=jj+1;
        end
        slope1=(deltaFTraces(ii,ii_pre+1-6+jj-1)-deltaFTraces(ii,ii_pre+1-6+jj-2))/dt_sample;
        taft=(jj-1)*dt_sample+((dFFo+0.5*(maxdFFo-dFFo)-deltaFTraces(ii,ii_pre+1-6+jj-2))/slope1);
        no_s_fwhm=no_s_fwhm+1;
        fwhm_per_spike(no_s_fwhm)=taft-tbef;
        
        %Plot the fwhm calculation
        
        try
            close 9
        catch
        end
        figure(9)
        
        plot([1:7+ii_post]*dt_sample,deltaFTraces(ii,ii_pre+1-6:end),'-k')
        hold on
        plot([1:10]*dt_sample,dFFo*ones(1,10),'-k')
        plot([maxii*dt_sample maxii*dt_sample],[min(deltaFTraces(ii,ii_pre+1-6:end)) max(deltaFTraces(ii,ii_pre+1-6:end))],'-r')
        plot([1:10]*dt_sample,(dFFo+0.5*(maxdFFo-dFFo))*ones(1,10),'-k')
        plot([tbef tbef],[min(deltaFTraces(ii,ii_pre+1-6:end)) max(deltaFTraces(ii,ii_pre+1-6:end))],'-b')
        plot([taft taft],[min(deltaFTraces(ii,ii_pre+1-6:end)) max(deltaFTraces(ii,ii_pre+1-6:end))],'-b')
    catch
    end
    pffft=1;
   
end

figure(8)
histogram(fwhm_per_spike,[0:1:15])

%Calculate roc
roc_data=[];

roc_data(1:no_resp,1)=dffPre;
roc_data(1:no_resp,2)=zeros(no_resp,1);

roc_data(1+no_resp:2*no_resp,1)=dffPeak;
roc_data(1+no_resp:2*no_resp,2)=ones(no_resp,1);

roc=roc_calc(roc_data);

save([inputPath(1:end-4) 'deltaFdivF.mat'],'deltaF_pre','deltaF_post','fwhm_per_spike','spike_time','spike_area_no')
pffft=1



