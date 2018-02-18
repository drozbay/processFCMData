function [ activeROIData ,bwActiveROI ] = getFCMActivity( imInput, inputROI, infoStruct, cfg )

% Get individual parameters
tWin = cfg.tWin;
pWin=tWin*4;
sdThresh = cfg.sdThresh;
minSize = cfg.minSize;
dt = infoStruct.dt;

xInput = size(imInput,2);
yInput = size(imInput,1);

numFrames = infoStruct.numFramesTrim;
frameStart = infoStruct.trimFrames(1);
frameEnd = infoStruct.trimFrames(2);
% Trim input image
imInputTrim = imInput(:,:,frameStart:frameEnd);

hWait = waitbar(0, sprintf('Filtering %d images...',numFrames));
for ii=1:numFrames
    imInputTrim(:,:,ii) = imgaussfilt(medfilt2(imInputTrim(:,:,ii),[1,1]),1);
    waitbar(ii/numFrames,hWait);
end
close(hWait);

%% Extract all areas from above
% Data structure for ROIs
fullROIData = struct(...
    'ROINum',[],...
    'intensity',zeros(numFrames,1),...
    'deltaF',zeros(numFrames,1),...
    'devVal',[],...
    'medianVal',[],...
    'hasEvent',[],...
    'numEvents',[],...
    'eventIdx',[],...
    'eventEdgeIdx',[],...
    'eventMaxIdx',[],...
    'deltaFPre',[],...
    'deltaFPeak',[],...
    'deltaFTraces',[]);

% Check if inputROI is labeled or not
if max(inputROI(:))>1
    % Labeled ROIs from input
    lwROI = inputROI;
    % Unlabeled ROIs
    bwROI = lwROI;
    bwROI(bwROI>0) = 1;
else
    % Unlabeled ROIs from input
    bwROI = logical(inputROI);
    % Labeled ROIs
    lwROI = bwlabel(inputROI,4);
end

% Get get mean intensity of each ROI for entire sequence
numROIs=length(regionprops(lwROI, imInputTrim(:,:,1), 'all'));
hWait = waitbar(0, sprintf('Processing %d images...',numFrames));
for ii = 1:numFrames
    tmpCurrentImage = imInputTrim(:,:,ii);
    ROIprop = regionprops(lwROI, tmpCurrentImage, 'all');
    for jj=1:numROIs
        fullROIData(jj).intensity(ii,1)=ROIprop(jj).MeanIntensity;
    end
    waitbar(ii/numFrames,hWait);
end
% Assign a number for each ROI
% Also get centroid locations
for jj = 1:numROIs
    fullROIData(jj).ROINum = jj;
    fullROIData(jj).Centroid = ROIprop(jj).Centroid;
end
clear tmp*;
close(hWait);

%% Processing of data
timeVector=linspace(0,(numFrames-1)*dt,numFrames);
% DeltaF/F Calculation
t0 = .5;
t1 = dt*3;
t2 = dt*6;
hWait = waitbar(0, sprintf('Performing deltaF/F on %d ROIs...',numROIs));
for jj = 1:numROIs
    fullROIData(jj).deltaF = deltaFCalc(timeVector,...
        fullROIData(jj).intensity,t0,t1,t2);
%     fullROIData(jj).deltaF = deltaFCalc(timeVector,...
%             15e-3+randi(2e4,length(fullROIData(jj).intensity),1)/1e6,t0,t1,t2);
    waitbar(jj/numROIs,hWait);
end
close(hWait);

% Get an estimate of the deviation
for jj=1:numROIs
    fullROIData(jj).devVal = mad(fullROIData(jj).deltaF,1);
    fullROIData(jj).medianVal = median(fullROIData(jj).deltaF);
end

% Find events using deviation threshold
isEvent=zeros(numFrames,numROIs);
for jj=1:numROIs
    devVal = fullROIData(jj).devVal;
    medVal = fullROIData(jj).medianVal;
    %Calculate events
    if devVal<0.1
        for ii=1:numFrames
            if fullROIData(jj).deltaF(ii) > ...
                    sdThresh*devVal+medVal
                isEvent(ii,jj)=1;
            end
        end
    end
    % Remove repeated vaues from isEvent vector
    tmpCleanEvents = [0; diff(isEvent(:,jj))];
    tmpCleanEvents(tmpCleanEvents<0) = 0;
    isEvent(:,jj) = isEvent(:,jj).*tmpCleanEvents;
    % Remove events that occur within tWin from a previous one
    for ii = 1:numFrames-tWin
       if isEvent(ii,jj)
           isEvent(ii+1:ii+tWin,jj) = 0;
       end
    end
end
clear tmp*;

fprintf('Total events: %d\n',numel(isEvent(isEvent==1)));

%% Count the number of events
for jj=1:numROIs
    if sum(isEvent(:,jj))>=1
        fullROIData(jj).hasEvent = 1;
        fullROIData(jj).eventIdx = find(isEvent(:,jj));
    else
        fullROIData(jj).hasEvent = 0;
    end
end
% % (OPTIONAL) Clear event with only one transient
% for jj = 1:numROIs
%     if length(fullROIData(jj).eventIdx)<2
%         fullROIData(jj).hasEvent = 0;
%     end    
% end
% Create new struct of ROIs with events
idxThr = [fullROIData.hasEvent]==1;
activeROIData = fullROIData(idxThr);
clear tmp*;

%% Get pre-transient and peak-transient values
for jj=1:length(activeROIData)
    kk = 0;
    numEvents = length(activeROIData(jj).eventIdx);
    deltaFPre = zeros(numEvents,1);
    deltaFPeaks = zeros(numEvents,1);
    deltaF = padarray(activeROIData(jj).deltaF,[tWin 0],'replicate');
    for ii=1:numEvents
        idxThr = activeROIData(jj).eventIdx(ii)+tWin;
        % Get the peak value tWin around the event
        deltaFPeaks(ii) = max(deltaF(idxThr-tWin:idxThr+tWin));
        % Find the index of the maximum value tWin around the event
        idxPeak = find(deltaF==deltaFPeaks(ii));
        if length(idxPeak) > 1
            warning('This should never happen');
            idxPeak = idxPeak(1);
        end
        % Find the index of the first value above 1/e*deltaFPeaks(ii)
        idxEdge = find(deltaF(idxThr-tWin:idxThr) >= ...
            deltaFPeaks(ii)*1/exp(1),1,'first')+idxThr-tWin-1;
        if isempty(idxEdge) || idxEdge>idxThr
            idxEdge = idxThr;
        end
        % Get the mean local value tWin ahead of the event edge
        deltaFPre(ii) = mean(deltaF(idxThr-tWin-1:idxThr-tWin+1));
        % Get windowed trace of transient, staying within edges
        if (idxPeak-tWin-pWin>0) && (idxPeak-tWin+pWin<numFrames)
            kk = kk + 1;
            deltaFTracesAtMax = deltaF(idxPeak-pWin:idxPeak+pWin);
            activeROIData(jj).deltaFTracesAtMax(:,kk) = deltaFTracesAtMax;
            deltaFTracesAtEdge = deltaF(idxEdge-pWin:idxEdge+pWin);
            activeROIData(jj).deltaFTracesAtEdge(:,kk) = deltaFTracesAtEdge;
            deltaFTraces = deltaF(idxThr-pWin:idxThr+pWin);
            activeROIData(jj).deltaFTraces(:,kk) = deltaFTraces;
        end
        % Update event index with index for the peak
        activeROIData(jj).eventEdgeIdx(ii) = idxEdge - tWin;
        activeROIData(jj).eventMaxIdx(ii) = idxPeak - tWin;
        fullROIData(activeROIData(jj).ROINum).eventEdgeIdx(ii) =...
            idxEdge - tWin;
        fullROIData(activeROIData(jj).ROINum).eventMaxIdx(ii) =...
            idxPeak - tWin;
    end
    activeROIData(jj).numEvents = length(deltaFPeaks);
    activeROIData(jj).deltaFPeak = deltaFPeaks;
    activeROIData(jj).deltaFPre = deltaFPre;
end

%% Display list of ROIs and indicate events
if ~exist('hfigAllTraces')
    hfigAllTraces = figure('units','inches','position',[5 0.5 9 9]);
elseif ~isgraphics(hfigAllTraces)
    hfigAllTraces = figure('units','inches','position',[5 0.5 9 9]);
end
figure(hfigAllTraces); clf;
subplot(1,2,1); cla;
imshow(bwROI.*0.5,'initialmagnification','fit'); hold on;
for jj = 1:length(fullROIData)
    idxThr = fullROIData(jj).ROINum;
    text(fullROIData(jj).Centroid(1),fullROIData(jj).Centroid(2),...
        sprintf('%d',idxThr),'Color','w');
    
end

subplot(1,2,2); cla;
hold on;
traceSep = 0.5;
for jj=1:numROIs
    devVal = fullROIData(jj).devVal;
    medVal = fullROIData(jj).medianVal;
    plot(timeVector, fullROIData(jj).deltaF - medVal+(jj-1)*traceSep);
    tmpEventIdx=ones(length(timeVector),1)*10e12;
    tmpEventIdx(fullROIData(jj).eventIdx) = 1;
    tmpEventMaxIdx=ones(length(timeVector),1)*10e12;
    tmpEventMaxIdx(fullROIData(jj).eventMaxIdx) = 1;
    tmpEventEdgeIdx=ones(length(timeVector),1)*10e12;
    tmpEventEdgeIdx(fullROIData(jj).eventEdgeIdx) = 1;
    plot(timeVector,tmpEventIdx/4+(jj-1)*traceSep,'.r','markersize',10);
    plot(timeVector,tmpEventMaxIdx/4+(jj-1)*traceSep,'.b','markersize',10);
    plot(timeVector,tmpEventEdgeIdx/4+(jj-1)*traceSep,'.g','markersize',10);
    plot(timeVector,repmat(devVal,numFrames,1)+(jj-1)*traceSep,...
        '--b','linewidth',0.1);
    plot(timeVector,repmat(sdThresh*devVal,numFrames,1)+(jj-1)*traceSep,...
        '--r','linewidth',0.1);
    plot(timeVector,zeros(1,numFrames)+(jj-1)*traceSep,...
        '-b','linewidth',0.1);
    text(-timeVector(end)/100,(jj-1)*traceSep,sprintf('%d',jj),...
        'HorizontalAlignment','Right');
end
hold off;
xlim([0 timeVector(end)]);
ylim([-traceSep (numROIs+1)*traceSep]);
set(gca,'ytick',[]);

%% Ask to remove ROIs manually
inputStr = input('Enter ROIs to remove (comma separated) [None]: ','s');
if isempty(inputStr)
    ROIsToRemove = [];
else
    c = textscan(inputStr,'%f','delimiter',',');
    ROIsToRemove = uint8(abs(round(c{1}')));
    disp('Will remove the following ROIs:')
    disp(ROIsToRemove);
end
%%
for jj=1:length(ROIsToRemove)
    curIdx = find([activeROIData.ROINum]==ROIsToRemove(jj),1,'first');
    activeROIData(curIdx) = [];
end
%% Save ROI file and data
% connComponents = bwconncomp(bwROI);
bwActiveROI = ismember(lwROI, [activeROIData.ROINum]);

[pname, fname, ~] = fileparts(infoStruct.fullPath);
% Write binary ROI file
imwrite(uint8(bwActiveROI),[pname,'/', fname, '_ROIEvents.tif'])
% Write colored ROI file
% roiCMap = hsv(length(activeROIData));
% colEventROIFilled = label2rgb(bwlabel(bwActiveROI), roiCMap, 'k');
% imwrite(colEventROIFilled,[pname,'/', fname, '_ROIEventsRGBFilled.png'])
% colEventROIEdges = label2rgb(bwlabel(bwperim(bwActiveROI)), roiCMap, 'k');
% imwrite(colEventROIEdges,[pname,'/', fname, '_ROIEventsRGBEdges.png'],...
%     'transparency',[0 0 0]);
% Save DeltaF Data file
save([pname,'/', fname, '_DeltaF.mat'],'',...
    'activeROIData','infoStruct','bwActiveROI');

end