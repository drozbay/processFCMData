%% processMouseFCM.m
% Baris Ozbay 1/13/2018
% Takes input of synchronized behavioral mouse video and mouse FCM imaging
% results (pre-processed to remove fiber artifact) and outputs information
% about the fluorescence transients and activity correlation
close all;
clear all;
%% User parameters
% Use config file for default parameters, or create it
if exist('processMouseFCM_cfg.mat', 'file') == 2
    load('processMouseFCM_cfg.mat');
else
    % Default cfg file parameters
    cfg.pnameFCMStart = pwd;
    cfg.pnameBehStart = pwd;
    cfg.useMOCO = 0;
    cfg.useROIFile = 0;
    cfg.useDeltaF = 0;
    cfg.startTime = 0;
    cfg.endTime = 0; % Note: Setting to 0 will default to end of series
    cfg.w = 20;
    cfg.tWin = 5;
    cfg.sdThresh = 6;
    cfg.minSize = 10;
    save('processMouseFCM_cfg.mat','cfg');
end

% Get the FCM image stack file
[fnameFCM,pnameFCM,nCancel] = uigetfile({[cfg.pnameFCMStart,'/*.tif']},...
    'Select the FCM file...');
if ~nCancel
    error('User cancelled');
end
% Get the mouse behavioral movie file
[fnameBeh,pnameBeh,nCancel] = uigetfile({[cfg.pnameBehStart,'/*.tif']},...
    'Select the Behavioral Movie file...');
if ~nCancel
    error('User cancelled');
end
% Get the motion correction file
% First check if there is a single csv file in directory
dirCsv = dir([pnameFCM,'/*.csv']);
if length(dirCsv) == 1
    fnameMOCO = dirCsv(1).name;
    pnameMOCO = dirCsv(1).folder;
    cfg.useMOCO = 1;
else
    % If there is not only one csv file in the FCM path, ask for file
    [fnameMOCO,pnameMOCO,nCancel] = uigetfile({[cfg.pnameFCMStart,'/*.csv']},...
        'Select the MOCO file (Optional)...');
    if ~nCancel
        cfg.useMOCO = 0;
    else
        cfg.useMOCO = 1;
    end
end

fullPathFCM = fullfile(pnameFCM,fnameFCM);
fullPathBeh = fullfile(pnameBeh,fnameBeh);
fullPathMOCO = [];
if cfg.useMOCO
    fullPathMOCO = fullfile(pnameMOCO,fnameMOCO);
end

% Save config file with path names
cfg.pnameFCMStart = pnameFCM;
cfg.pnameBehStart = pnameBeh;
save('processMouseFCM_cfg.mat','cfg');

% If *ROI.tif file exists, request to use it or generate new one 
[~,nameFCM,~] = fileparts(fullPathFCM);
if exist(fullfile(pnameFCM,[nameFCM,'_ROI.tif']),'file')==2
    tmpInputStr = input('Use existing *ROI.tif file? Y/N [Y]: ','s');
    if isempty(tmpInputStr)
        tmpInputStr = 'Y';
    end
    tmpInputStr = lower(tmpInputStr);
    if strcmp(tmpInputStr(1),'n')
        cfg.useROIFile = 0;
    else
        cfg.useROIFile = 1;
    end
    clear tmp*;
else
    cfg.useROIFile = 0;
end

% If *DeltaF.mat file exists, request to use it or generate new one 
[~,nameFCM,~] = fileparts(fullPathFCM);
if (exist(fullfile(pnameFCM,[nameFCM,'_DeltaF.mat']),'file')==2)...
    && cfg.useROIFile
    % If DeltaF.mat file exists and user has chosen to use ROI.tif file
    tmpInputStr = input('Use existing *DeltaF.mat file? Y/N [Y]: ','s');
    if isempty(tmpInputStr)
        tmpInputStr = 'Y';
    end
    tmpInputStr = lower(tmpInputStr);
    if strcmp(tmpInputStr(1),'n')
        cfg.useDeltaF = 0;
    else
        cfg.useDeltaF = 1;
    end
    clear tmp*;
else
    cfg.useDeltaF = 0;
end
save('processMouseFCM_cfg.mat','cfg');

%% Get info from FCM and Behavioral Movie files
% FCM file:
infoFCM = getTimeSeriesInfo(fullPathFCM); 
% Behavioral movie file info:
infoBeh = getTimeSeriesInfo(fullPathBeh);
% Rescale framerates so the durations match
if abs(infoFCM.duration/infoBeh.duration-1)>0.02
    warning('Difference in acquisition durations is larger than 2%');
    warning('FCM file duration: %.1f sec',infoFCM.duration);
    warning('Beh file duration: %.1f sec',infoBeh.duration);
end
infoBeh.dt = infoBeh.dt*infoFCM.duration/infoBeh.duration;
infoBeh.duration = infoBeh.dt*(infoBeh.numFrames-1);

% If default end time selected, assign it to duration of acquisition
if cfg.endTime == 0
    cfg.endTime = infoBeh.duration;
end

% Set trimFrames parameters based on configured range
% Find matching frames in FCM stack
infoBeh.trimFrames = round(...
    ([cfg.startTime cfg.endTime]+infoBeh.dt)/infoBeh.dt );
% Find matching frames in FCM stack
infoFCM.trimFrames = round(...
    ([cfg.startTime cfg.endTime]+infoFCM.dt)/infoFCM.dt );

%% Import previously saved info for the time range
if cfg.useDeltaF
    % If user selected to use DeltaF file then use that time range
    % Load DeltaF file from folder
    load([pnameFCM,'/',nameFCM,'_DeltaF.mat']);
    % Change existing frames to match external file if necessary
    if max(infoStruct.trimFrames ~= infoFCM.trimFrames)
        infoFCM.trimFrames = infoStruct.trimFrames;
        cfg.startTime = (infoStruct.trimFrames(1)-1)*infoFCM.dt;
        cfg.endTime = (infoStruct.trimFrames(2)-1)*infoFCM.dt;
        fprintf(['Saved time range does not match configured range',...
            newline,'Using range %.1f to %.1f sec from existing file',...
            newline],...
            cfg.startTime,...
            cfg.endTime);
    end
    % Remove infoStruct (only used to get frame range)
    clear infoStruct;
    changeTrim = 0; % Don't request new range from user
else
    % If not using DeltaF fill ask to change the
    % existing time window for analysis
    fprintf(['Currently configured time window:',newline,...
        'startTime: %.1f sec',newline,...
        'endTime: %.1f sec',newline],...
        cfg.startTime,cfg.endTime);
    tmpInputStr = input('Use this window? Y/N [Y]: ','s');
    if isempty(tmpInputStr)
        tmpInputStr = 'Y';
    end
    tmpInputStr = lower(tmpInputStr);
    if strcmp(tmpInputStr(1),'y')
        changeTrim = 0;
    else
        changeTrim = 1;
    end
    clear tmp*;
end

% If start or end times don't make sense, require user to input new numbers
if cfg.startTime>infoBeh.duration || cfg.endTime>infoBeh.duration
    disp('Selected time range exceeds duration of imaging');
    changeTrim = 1;
end

save('processMouseFCM_cfg.mat','cfg');

%% Import image frames
imFCM = importImageSeries(infoFCM);
imBeh = importImageSeries(infoBeh);

%% Trim in selected time range 
if changeTrim
    % If user needs to select trimming parameters, open a dialog box with
    % behavioral movie stack to select start and end frames
    infoBeh = getTrimFrames(imBeh,infoBeh);
    cfg.startTime = (infoBeh.trimFrames(1)-1)*infoBeh.dt;
    cfg.endTime = (infoBeh.trimFrames(2)-1)*infoBeh.dt;
end

% Find matching frames in FCM stack
infoBeh.trimFrames = round(...
    ([cfg.startTime cfg.endTime]+infoBeh.dt)/infoBeh.dt );
% Find matching frames in FCM stack
infoFCM.trimFrames = round(...
    ([cfg.startTime cfg.endTime]+infoFCM.dt)/infoFCM.dt );

% Make sure frames fall in range limits
if infoBeh.trimFrames(1) < 1
    infoBeh.trimFrames(1) = 1;
end
if infoBeh.trimFrames(2) > infoBeh.numFrames
    infoBeh.trimFrames(2) = infoBeh.numFrames;
end
if infoFCM.trimFrames(1) < 1
    infoFCM.trimFrames(1) = 1;
end
if infoFCM.trimFrames(2) > infoFCM.numFrames
    infoFCM.trimFrames(2) = infoFCM.numFrames;
end
% Trim down images and update info structs
infoBeh.numFramesTrim = infoBeh.trimFrames(2)-infoBeh.trimFrames(1)+1;
infoBeh.durationTrim = infoBeh.dt*(infoBeh.numFramesTrim-1);
infoFCM.numFramesTrim = infoFCM.trimFrames(2)-infoFCM.trimFrames(1)+1;
infoFCM.durationTrim = infoFCM.dt*(infoFCM.numFramesTrim-1);

save('processMouseFCM_cfg.mat','cfg');

%% Set other parameters
% Calculated parameters:
% Number of pixels for square window to determine event areas
cfg.w = roundOdd(mean(size(imFCM,1),size(imFCM,2))/10);
% Minimum region size
cfg.minSize = round(mean(size(imFCM,1),size(imFCM,2))/10);

% Defaults for configurable parameters
if ~isfield(cfg,'tWin') || ~isfield(cfg,'sdThresh')
    % If any field is missing, revert to defaults
    cfg.tWin = 5; % Time window
    cfg.sdThresh = 5; % SD threshold for event
end

% Show dialog box so user may choose different parameters
prompt = {'Enter minimum ROI area:',...
    'Enter time window (tWin):',...
    'Enter signal threshold (sdThresh):'};
num_lines = 1;
tmpDefaultAns = {num2str(cfg.minSize),num2str(cfg.tWin),...
    num2str(cfg.sdThresh)};
tmpAnswer = inputdlg(prompt,'Input parameters',1,tmpDefaultAns);
if isempty(tmpAnswer)
    error('User cancelled');
end
cfg.minSize = str2double(tmpAnswer{1});
cfg.tWin = str2double(tmpAnswer{2});
cfg.sdThresh = str2double(tmpAnswer{3});
clear tmp*;

% Update cfg
save('processMouseFCM_cfg.mat','cfg');

%% Get ROI info from FCM file
if cfg.useROIFile
    % If existing ROI file is used, read it in
    bwROI=logical(imread([pnameFCM,'/',nameFCM,'_ROIManual.tif']));
    [bwROI, wsROI, lwROI, rgbROI] = getFCMROI(imFCM,infoFCM,cfg, bwROI);
else
    % If ROI file is not used, run getFCMROI
    [bwROI, wsROI, lwROI, rgbROI] = getFCMROI(imFCM,infoFCM,cfg);
end

% Plot ROIs
hfigROI = figure(14); clf;
set(hfigROI,'units','inches','position',[4,4,8,3])
subplot(1,3,1);
imshow(bwROI);
title('Active area mask');
% Watershed
subplot(1,3,2);
imshow((double(wsROI)+double(bwROI))/2);
title('Segmented ROIs');
% Color plot
subplot(1,3,3);
imshow(rgbROI);
title('Individual ROIs');

%% Process ROIs and recover activity of regions
if ~cfg.useDeltaF
    % If DeltaF file is not used, run getFCMActivity
    [activeROIData, bwActiveROI] = getFCMActivity( imFCM, bwROI, infoFCM, cfg );
end

%% Get optical movement data from behavioral movie
behSubFactor = 5;
[~,flowMag1D] = getBehOpticalFlow(imBeh,infoBeh,behSubFactor);
behData = movmedian(flowMag1D,50);
behTime = linspace(0,(infoFCM.numFramesTrim-1)*infoFCM.dt,...
    length(behData));
figure(101); plot(behData);

%% Prepare required quantities for plotting
% Make ROI numbers sequential
for jj = 1:length(activeROIData)
    activeROIData(jj).ROINum = jj;
end
ROINum = horzcat(activeROIData.ROINum);

dfLow = 0.3;
numFrames = infoFCM.numFramesTrim;
deltaF = horzcat(activeROIData.deltaF);
% Array of highest deltaF event value for each ROI
deltaFPeakROI = zeros(1,length(activeROIData));
% Array of number of events for each ROI
deltaFNumEvent = zeros(1,length(activeROIData));
for jj = 1:length(deltaFPeakROI)
    deltaFPeakROI(jj) = max(activeROIData(jj).deltaFPeak);
    deltaFNumEvent(jj) = length(activeROIData(jj).deltaFPeak);
end
% Quantities based on each individual event:
deltaFTraces = horzcat(activeROIData.deltaFTracesAtEdge);
eventTimes = horzcat(activeROIData.eventEdgeIdx)*infoFCM.dt;
deltaFPeak = vertcat(activeROIData.deltaFPeak);
% Array of numbers corresponding to ROI of each event
eventROINums = [];
for ii = 1:length(activeROIData)
    for jj=1:length(activeROIData(ii).eventIdx)
        eventROINums = [eventROINums activeROIData(ii).ROINum];
    end
end
pWin = cfg.tWin*4;
% Time axis
timeVector=linspace(0,(numFrames-1)*infoFCM.dt,numFrames);

% Colormap
roiCMap = hsv(round(length(activeROIData)*1.5));
roiCMap = roiCMap(1:length(activeROIData),:)*0.9;

% Grid colors
gridCol = [0.7 0.7 0.7];

% ROI images
bwlEvent = bwlabel(bwActiveROI);
bwEventRGB = label2rgb(bwlEvent, roiCMap, 'k');
bwlEventEdgesRGB = label2rgb(bwlabel(bwperim(bwlEvent)), roiCMap, 'k');

% Histogram edges
histTimeWindow = 20;  % seconds
histEdgesFCM = (0:histTimeWindow:timeVector(end)+histTimeWindow);

% Directory to save outputs
[pnameOut, fnameOut, ~] = fileparts(infoFCM.fullPath);
[~,~,~] = mkdir([pnameFCM,'/Outputs']);

% Output figure properties
fontsize = 6;
dpiOut = '-r600';

%% Plot summary
hfigSumm = figure(3); clf;
set(hfigSumm,'units','inches','position',[2,2,10,8]);
% RGB ROI image
subplot(2,2,1); cla;
imshow(bwEventRGB*0.8,'initialmagnification','fit'); hold on;
for jj = 1:length(activeROIData)
    tmpIdx = activeROIData(jj).ROINum;
    text(activeROIData(jj).Centroid(1),activeROIData(jj).Centroid(2),...
        sprintf('%d',tmpIdx),'Color','w');
end

subplot(2,2,2); cla; hold on;
tmpTraceSep = 0.5;
for jj=1:length(activeROIData)
    tmpStdVal = activeROIData(jj).devVal;
    tmpMedVal = activeROIData(jj).medianVal;
    plot(timeVector,(activeROIData(jj).deltaF-tmpMedVal) +...
        (jj-1)*tmpTraceSep,...
        '-','Color',roiCMap(jj,:));
    text(-timeVector(end)/100,(jj-1)*tmpTraceSep,...
        sprintf('%d',activeROIData(jj).ROINum),...
        'Color',roiCMap(jj,:),'HorizontalAlignment','Right');
    % Add arrows at events
    for ii=1:length(activeROIData(jj).eventEdgeIdx)
        tmpTimePoint = activeROIData(jj).eventEdgeIdx(ii)*infoFCM.dt;
        tmpH = annotation('arrow');
        set(tmpH,'parent', gca, ...
            'position', [tmpTimePoint (jj-1.5)*tmpTraceSep 0 tmpTraceSep*0.5], ...
            'HeadLength', 2, 'HeadWidth', 2, 'HeadStyle', 'plain', ...
            'Color', roiCMap(jj,:)*0.8, 'LineWidth', 0.5);
    end
end
xlim([0 timeVector(end)]);
ylim([-tmpTraceSep (length(activeROIData)+1)*tmpTraceSep]);
set(gca,'YTick',[]);
xlabel('Time (sec)');
ylabel(['Region #',newline]);
set(gca,'XGrid','on');
set(gca,'XMinorGrid','on');

% Event windowed timecourses
subplot(2,2,3);
hold on;
for ii = 1:size(deltaFTraces,2)
    plot((0:pWin*2)*infoFCM.dt,deltaFTraces(:,ii),'-','color',roiCMap(eventROINums(ii),:),'linewidth',0.5)
end
plot([pWin+1 pWin+1]*infoFCM.dt,[min(deltaFTraces(:)) max(deltaFTraces(:))],'-r')
plot((0:pWin*2)*infoFCM.dt,mean(deltaFTraces,2),'-k','linewidth',2)

xlabel('Time (seconds)');
title('Transient event timecourses');
xlim([0 pWin*2]*infoFCM.dt);
ylabel('DeltaF/F');
set(gca,'XGrid','on');
set(gca,'XMinorGrid','on');

% Show spikes raster plot
subplot(2,2,4); cla; colormap(roiCMap);
scatter(eventTimes,deltaFPeak,50,eventROINums,'s','filled');
xlim([0,timeVector(end)]);
title('Calcium Transients');
xlabel('Time (sec)');
ylabel('DeltaF/F');
set(gca,'XGrid','on');
set(gca,'XMinorGrid','on');
% colorbar;
clear tmp*;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% % Make output figures % %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sc = 1; % Plot scale
%% Scatter plot
dfLim = 0;
limIdx = deltaFPeak>dfLim;
deltaFPeaksLim = deltaFPeak(limIdx);
eventTimesLim = eventTimes(limIdx);

hfigSct = figure(11); clf;
set(hfigSct,'Color',[1 1 1]);
set(hfigSct,'units','inches','position',[3,9,1.25*sc,0.3*sc]);
colormap(roiCMap);

% % Scatter plot when y axis is deltaF
% sct = scatter(eventTimesLim,deltaFPeaksLim,5*sc,eventROINums,...
%     's','filled');
% set(gca,'color',[1 1 1]);
% set(gca,'gridcolor',gridCol,'gridalpha',1);
% set(gca,'minorgridcolor',gridCol,'minorgridalpha',0.25,...
%     'MinorGridLineStyle','-');
% xlim([0,timeVector(end)]);
% ylim([0,1.05]);
% ylim([0,max(eventROINums)+2]);
% box on;
% ylabel('ROI Color');
% set(gca,'XTick',[0:100:1000]);
% set(gca,'Ytick',[]);
% set(gca,'XGrid','on');
% set(gca,'XMinorGrid','on');
% set(gca,'Units','Inches','Position',[0.005 0.005 1.24*sc 0.29*sc]);
% set(gca,'FontSize',fontsize);
% [~,~,~] = mkdir([pnameFCM,'/Outputs']);
% print([pnameOut,'/Outputs/',fnameOut,'_ScatterPlot'],...
%     '-dtiff', dpiOut);

% Scatter plot when y axis is roiNum
sct = scatter(eventTimesLim,eventROINums,1,eventROINums,...
    's','filled');
hold on;
scatter(eventTimesLim+1,eventROINums,1,eventROINums,...
    's','filled');
scatter(eventTimesLim+4,eventROINums,1,eventROINums,...
    's','filled');
scatter(eventTimesLim+7,eventROINums,1,eventROINums,...
    's','filled');

set(gca,'color',[1 1 1]);
set(gca,'gridcolor',gridCol,'gridalpha',1);
set(gca,'minorgridcolor',gridCol,'minorgridalpha',0.25,...
    'MinorGridLineStyle','-');
xlim([0,timeVector(end)]);
ylim([0,max(eventROINums)+1]);
box on;
set(gca,'Ytick',[]);
set(gca,'XGrid','on');
set(gca,'XMinorGrid','on');
set(gca,'Units','Inches','Position',[0.005 0.005 1.24*sc 0.29*sc]);
set(gca,'FontSize',fontsize);
[~,~,~] = mkdir([pnameFCM,'/Outputs']);
print([pnameOut,'/Outputs/',fnameOut,'_ScatterPlot1'],...
    '-dtiff', dpiOut);


% Scatter plot when y axis is ROI number
% colormap('lines');
% sct = scatter(eventTimesLim,eventROINums,5,eventROINums,...
%     'o','filled');
% set(gca,'color',[1 1 1]);
% set(gca,'gridcolor',gridCol,'gridalpha',1);
% set(gca,'minorgridcolor',gridCol,'minorgridalpha',0.25,...
%     'MinorGridLineStyle','-');
% xlim([0,timeVector(end)]);
% ylim([0,max(eventROINums)+2]);
% box on;
% set(gca,'Ytick',[]);
% set(gca,'XGrid','on');
% set(gca,'XMinorGrid','on');
% set(gca,'Units','Inches','Position',[0.005 0.005 1.24*sc 0.29*sc]);
% set(gca,'FontSize',fontsize);
% [~,~,~] = mkdir([pnameFCM,'/Outputs']);
% print([pnameOut,'/Outputs/',fnameOut,'_ScatterPlot2'],...
%     '-dtiff', dpiOut);
% 
%% Events histogram
hfigHist = figure(12); clf;
set(hfigHist,'Color',[1 1 1]);
set(hfigHist,'units','inches','position',[3,7.5,1.25*sc,0.3*sc]);

hHist = histogram([activeROIData.eventEdgeIdx]*infoFCM.dt,histEdgesFCM);
maxCounts = max(hHist.BinCounts);
set(gca,'color',[1 1 1]);
set(gca,'gridcolor',gridCol,'gridalpha',1);
set(gca,'minorgridcolor',gridCol,'minorgridalpha',0.25,...
    'MinorGridLineStyle','-');
xlim([0,timeVector(end)]);
ylim([0,maxCounts+1]);
box on;
set(gca,'XTick',[0:100:1000]);

set(gca,'Ytick',[0 maxCounts]);
set(gca,'XGrid','on');
set(gca,'XMinorGrid','on');
set(gca,'TickLength',[0 0]);
set(gca,'Units','Inches','Position',[0.005 0.005 1.24*sc 0.29*sc]);
% xlabel('Time (sec)');
% ylabel('DeltaF/F');
% colorbar;
set(gca,'FontSize',fontsize);
% 
print([pnameOut,'/Outputs/',fnameOut,'_EventHist'],...
   '-dtiff', dpiOut);

%% Behavioral data plot
% Behavioral data
behMean = circshift(movmean(behData,histTimeWindow/infoBeh.dt/10),0,1);

hfigBeh = figure(13); clf;
set(hfigBeh,'Color',[1 1 1]);
set(hfigBeh,'units','inches','position',[3,6,1.25*sc,0.3*sc]);
% behMean = movmax(abs(gradient(behData)),20);
plot(behTime,behMean,'-k','linewidth',1);

set(gca,'color',[1 1 1]);
set(gca,'gridcolor',gridCol,'gridalpha',1);
set(gca,'minorgridcolor',gridCol,'minorgridalpha',0.25,...
    'MinorGridLineStyle','-');
xlim([0,timeVector(end)]);
ylim([0,max(behMean)*1.1]);
box on;
set(gca,'XTick',[0:100:1000]);
set(gca,'Ytick',[]);
set(gca,'XGrid','on');
set(gca,'XMinorGrid','on');
set(gca,'Units','Inches','Position',[0.005 0.005 1.24*sc 0.29*sc]);
% xlabel('Time (sec)');
% ylabel('DeltaF/F');
% colorbar;
set(gca,'FontSize',fontsize);

print([pnameOut,'/Outputs/',fnameOut,'_BehavioralPlot'],...
   '-dtiff', dpiOut);

%% Subset of traces plot
numToShow = 2;
% Sort by a weight of both events and peak:
weight = 0.5;
[~,deltaFSortIdx] = sort(...
    weight*deltaFNumEvent./max(deltaFNumEvent) + ...
    (1-weight)*deltaFPeakROI./max(deltaFPeakROI),...
    'descend');

deltaFSorted = deltaF(:,deltaFSortIdx);
ROINumSorted = ROINum(:,deltaFSortIdx);
ROINumSubset = sort(ROINumSorted(:,1:numToShow));

% Choose ROIs manually
% ROINumSubset = [ROINumSubset 3 11 22]; %For 13
% ROINumSubset = [3,4,9,15,24]; %For 17
% ROINumSubset = [1,4,8,17,36]; %For 11

ROINumSubset = sort(ROINumSubset);
 
deltaFSubset = deltaF(:,ROINumSubset);
hfigTraces = figure(14); clf;
set(hfigTraces,'Color',[1 1 1]);
set(hfigTraces,'units','inches','position',[3,0.5,3.3*sc,1*sc]);

hold on;
traceSep = 0.5; 
for jj = 1:length(ROINumSubset)
    curROIIdx = ROINumSubset(jj);
    plot(timeVector,deltaFSubset(:,jj)*2+(jj-1)*traceSep,...
        '-','color',roiCMap(ROINumSubset(jj),:),'linewidth',0.25);
    % Place current line plot at the top
    chH = get(gca,'Children');
    set(gca,'Children',[chH(2:end);chH(1)])
    % Add numbers for ROIs
    text(-5,(jj-1)*traceSep,...
        sprintf('%d',activeROIData(jj).ROINum),'fontsize',fontsize*0.8,...
        'Color',roiCMap(curROIIdx,:),'HorizontalAlignment','Right');
%     text(0,(jj-1)*traceSep,...
%         sprintf('%d',curROIIdx),'fontsize',fontsize,...
%         'Color',roiCMap(curROIIdx,:),'HorizontalAlignment','Right');
    for ii=1:length(activeROIData(curROIIdx).eventEdgeIdx)
        tmpTimePoint = activeROIData(curROIIdx).eventEdgeIdx(ii)*infoFCM.dt;
        tmpH = annotation('arrow');
        set(tmpH,'parent', gca, ...
            'position', [tmpTimePoint (jj-1.3)*traceSep 0 traceSep*0.1], ...
            'HeadLength', 1.5, 'HeadWidth', 1.5, 'HeadStyle', 'plain', ...
            'Color', roiCMap(curROIIdx,:)*0.8, 'LineWidth', 0.5);
    end
end
set(gca,'color',[1 1 1]);
% xlim([0 pWin*2+1]*infoFCM.dt);
xlim([0 timeVector(end)]);
ylim([-traceSep*1 traceSep*(length(ROINumSubset))+0.5]);
box off;
% set(gca,'Xtick',[0 5 10 15]);
set(gca,'TickDir','out');
set(gca,'YColor',[1 1 1]);
% set(gca,'XGrid','on');
% set(gca,'XMinorGrid','on'); 
% curPos = get(gca,'Position');
set(gca,'Units','Inches','Position',[0.2 0.05 3.09*sc 0.85*sc]);
set(gca,'FontSize',fontsize);

print([pnameOut,'/Outputs/',fnameOut,'_deltaFPlot'],...
    '-dtiff', '-r900');

%% Mean traces plot
hfigWind = figure(15); clf;
set(hfigWind,'Color',[1 1 1]);
set(hfigWind,'units','inches','position',[9,7,1*sc,1.5*sc]);
hold on;
tracesTimeAxis = ((0:pWin*2)-pWin*0.5)*infoFCM.dt;
for ii = 1:size(deltaFTraces,2)
    plot(tracesTimeAxis,deltaFTraces(:,ii),'-',...
        'color',[0.8 0.8 0.8],'linewidth',0.5)
end
    plot(tracesTimeAxis,mean(deltaFTraces,2),'-k','linewidth',2)
set(gca,'color',[1 1 1]);
xlim([0 20]);
ylim([-.2,1]);
box off;
set(gca,'Xtick',[0 5 10 15 20 25 30 40]);
set(gca,'Ytick',[0 0.5 1 1.5]);
set(gca,'TickDir','out');
% set(gca,'YColor',[1 1 1]);
% set(gca,'XGrid','on');
% set(gca,'XMinorGrid','on');
% set(gca,'Units','Inches','Position',[0.2 0.18 0.75 0.75]);
set(gca,'FontSize',fontsize);
 
print([pnameOut,'/Outputs/',fnameOut,'_TracesPlot'],...
    '-dtiff', dpiOut);


%% ROI plots
% Without numbers
hfigROI = figure(16); clf;
set(hfigROI,'Color',[1 1 1]);
set(hfigROI,'units','inches','position',[9,4,1*sc,1*sc]);
bwlEventEdgesRGB = label2rgb(bwlabel(bwperim(bwlEvent)), lines(), 'k');

% With all regions labeled
imshow(bwlEventEdgesRGB,'initialmagnification','fit'); hold on;
for jj = 1:length(activeROIData)
    tmpIdx = activeROIData(jj).ROINum;
    text(activeROIData(jj).Centroid(1)-5,activeROIData(jj).Centroid(2),...
        sprintf('%d',tmpIdx),'Color','w','fontsize',6);
end
print([pnameOut,'/Outputs/',fnameOut,'_rgbROI_allLabel'],...
    '-dtiff', dpiOut);

% With subset of regions labeled corresponding to traces shown
imshow(bwlEventEdgesRGB,'initialmagnification','fit'); hold on;
for jj = 1:length(ROINumSubset)
    tmpIdx = ROINumSubset(jj);
    text(activeROIData(tmpIdx).Centroid(1),...
        activeROIData(tmpIdx).Centroid(2),...
        sprintf('%d',jj),'Color','w','FontSize',fontsize);
end
print([pnameOut,'/Outputs/',fnameOut,'_rgbROI_tracesLabel'],...
    '-dtiff', dpiOut);
imwrite(bwlEventEdgesRGB,[pnameOut,'/Outputs/',fnameOut,...
    '_rgbROI_noLabel.png'],'transparency',[0 0 0]);

%% Motion regression
% Size of window for binning
edgeWindow = 20;  % seconds
% Lower limit of peak deltaF/F events to use
lowerLim = 0;
% Create arrays for binning
edgesForBehBinning = (0:edgeWindow:timeVector(end)+edgeWindow);
eventTimesBinned = histcounts(eventTimesLim,edgesForBehBinning);
behBinned = zeros(1,length(eventTimesBinned));
deltaFPeaksBinned = zeros(1,length(eventTimesBinned));
% Average out behavioral data for bin sampling
behMean = circshift(movmax(behData,histTimeWindow/infoBeh.dt/10),0,1)/5-5;
% Sample behavioral and deltaF data for each bin window
kk = 0;
for ii = 1:length(eventTimesBinned)
   curTime = edgesForBehBinning(ii)+edgeWindow/2;
   behIdx = find(behTime >= curTime,1);
   if isempty(behIdx)
      behIdx = length(behMean);
   end
   behBinned(ii) = behMean(behIdx);
   % Filter out data that is below the deltaF/F limit
   deltaFIdx = (eventTimesLim>=edgesForBehBinning(ii))...
       & (eventTimesLim<=edgesForBehBinning(ii)+edgeWindow);
   deltaFPeaksCur = deltaFPeaksLim(deltaFIdx);
   for jj = 1:length(deltaFPeaksCur)
       if deltaFPeaksCur(jj)<lowerLim
            deltaFPeaksCur(jj) = 0;
       end
   end
   deltaFPeaksCur(~deltaFPeaksCur) = [];
   if ~isempty(deltaFPeaksCur)
       deltaFPeaksBinned(ii) = max(deltaFPeaksCur);   
   else
       deltaFPeaksBinned(ii) = 0;
   end
end
% % Remove extra values from arrays
behBinned(~deltaFPeaksBinned) = [];
eventTimesBinned(~deltaFPeaksBinned) = [];
deltaFPeaksBinned(~deltaFPeaksBinned) = [];

% Linear Regression
Y = eventTimesBinned';
% Y = deltaFPeaksBinned';
X = behBinned';
XX = [ones(length(X),1) X];
b = XX\Y;
regY = XX*b;

[rhoBehCorr,pBehCorr] = corr(X,Y);

% Without numbers
hfigBehReg = figure(171); clf;
set(hfigBehReg,'Color',[1 1 1]);
set(hfigBehReg,'units','inches','position',[11,1,3,2]);
colormap(jet());

hold on;
scatter(X,Y,40,deltaFPeaksBinned,'s','filled');
% scatter(behBinned,deltaFPeaksBinned,30,eventTimesBinned,'o','filled');
maxY = max(eventTimesBinned);
minY = min(eventTimesBinned);
plot(X,regY,'-k','linewidth',2);
set(gca,'TickDir','out',...
    'XTick',[0,max(behBinned)/2,max(behBinned)],...
    'YTick',[0,maxY/2,maxY]);
xtickformat('%d')
xlim([0,max(behBinned)]);
ylim([0,maxY]);
ylabel('Number of Events');
xlabel('Motion (mm/s)');
text(max(behBinned)*0.05,...
    maxY*1.1,...
    sprintf('\\rho = %.3f\np_{sign} = %.3f',rhoBehCorr,pBehCorr),...
    'HorizontalAlignment','left','VerticalAlignment','top',...
    'FontSize',10);
hc = colorbar;
set(get(hc,'ylabel'),'string','\DeltaF/F');

print([pnameOut,'/Outputs/',fnameOut,'_BehavioralScatter'],...
    '-dtiff', dpiOut);

%% Plot histograms of thresholded deltaF
figure(172); clf; colormap('jet'); hold on;
thresh = 0.3;
behEdges = linspace(0,max(behBinned),20);

h1 = histogram(behBinned(deltaFPeaksBinned<thresh),behEdges);
h2 = histogram(behBinned(deltaFPeaksBinned>thresh),behEdges);
% h2.BinWidth = h1.BinWidth;
xlabel('Motion (mm/s)');
ylabel('Number of events');
legend(sprintf('\\DeltaF/F < %.2f',thresh),...
    sprintf('\\DeltaF/F > %.2f',thresh));

%% Motion correction
if cfg.useMOCO
    mocoData = importMOCO(fullPathMOCO,infoFCM);
end
hfigMOCO = figure(18); clf;
set(hfigMOCO,'Color',[1 1 1]);
set(hfigMOCO,'units','inches','position',[3,3,1.01*sc,0.3*sc]);
plot(mocoData.t_s,mocoData.mag_um);

set(gca,'color',[1 1 1]);
set(gca,'gridcolor',gridCol,'gridalpha',1);
set(gca,'minorgridcolor',gridCol,'minorgridalpha',0.25,...
    'MinorGridLineStyle','-');
xlim([0,timeVector(end)]);
% ylim([0,max(behMean)*1.1]);
box on;
set(gca,'XTick',[0:100:1000]);
set(gca,'XGrid','on');
set(gca,'XMinorGrid','on');
set(gca,'fontsize',fontsize);


set(gca,'Units','Inches','Position',[0.2 0.05 1*sc 0.2*sc]);
print([pnameOut,'/Outputs/',fnameOut,'_MOCOPlotNumbers'],...
   '-dtiff', dpiOut);

set(gca,'Ytick',[]);
set(gca,'Units','Inches','Position',[0.005 0.005 1.24*sc 0.29*sc]);
print([pnameOut,'/Outputs/',fnameOut,'_MOCOPlot'],...
   '-dtiff', dpiOut);
%% Get the max projection images of time points with events 
rWin = 1; % Odd numbered window around peak event to capture frames
[imFCM_maxProjActive, imFCM_wEvents, imFCM_noEvents] = ...
    getActiveFrames( imFCM, activeROIData, infoFCM, rWin );

hfigMaxProj = figure(19); clf;
set(hfigMaxProj,'Units','Inches','Position',[12 5 3 3]); 
colormap('hot');
imagesc(imFCM_maxProjActive);
title('Max projection');

% Save Max Projection file
imwrite(uint16(imFCM_maxProjActive./max(imFCM_maxProjActive(:)).*2^16-1),...
    [pnameOut,'/Outputs/', fnameOut, '_maxProj.tif']);

% Write images with separated events frames
maxSignal = max(imFCM_wEvents(:));
imFCM_wEvents = uint16(imFCM_wEvents./maxSignal.*2^16-1);
imwrite(imFCM_wEvents(:,:,1),[pnameOut,'/', fnameOut, '_wEvents.tif']);
for ii = 2:size(imFCM_wEvents,3)
    imwrite(imFCM_wEvents(:,:,ii),[pnameOut,'/', fnameOut, '_wEvents.tif'],...
        'writemode','append');
end

imFCM_noEvents = uint16(imFCM_noEvents./maxSignal.*2^16-1);
imwrite(imFCM_noEvents(:,:,1),[pnameOut,'/', fnameOut, '_noEvents.tif']);
for ii = 2:size(imFCM_noEvents,3)
    imwrite(imFCM_noEvents(:,:,ii),[pnameOut,'/', fnameOut, '_noEvents.tif'],...
        'writemode','append');
end

%% Get the fwhm
numEvents = size(deltaFTraces,2);
winLength = size(deltaFTraces,1);
winStart = floor(winLength/2.5);
winEnd = ceil(winLength*0.7);
eventFWHM=zeros(1,numEvents);
deltaFPre=zeros(1,numEvents);
deltaFMax=zeros(1,numEvents);
idxBefInt=zeros(1,numEvents);
idxAftInt=zeros(1,numEvents);
for ii=1:numEvents
    deltaFPre(ii) = deltaFTraces(winStart,ii);

    %Get the baseline
    dFBase=0;%median(deltaFTraces(1:winStart,ii));

    %Find the maximum
    [deltaFMax(ii), maxIdx]=max(deltaFTraces(winStart:winEnd,ii));
    maxIdx = maxIdx + winStart - 1;
  
    %Find points around the half way before and after marks
    halfMax = dFBase+0.5*(deltaFMax(ii)-dFBase);
    idxBef = find(deltaFTraces(winStart:winEnd,ii) > halfMax,1,'first')+winStart-2;
    idxAft = find(deltaFTraces(winStart:winEnd,ii) > halfMax,1,'last')+winStart;
    
    % Find the slope around each point
    slopeBef = deltaFTraces(idxBef,ii) - deltaFTraces(idxBef-1,ii);
    slopeAft = deltaFTraces(idxAft,ii) - deltaFTraces(idxAft-1,ii);
    % Interpolate the time point of the actual half max
    x0 = idxBef;
    x1 = idxBef+1;
    y0 = deltaFTraces(idxBef,ii);
    y1 = deltaFTraces(idxBef+1,ii);
    idxBefInt(ii) = (x0*(y1-halfMax)+x1*(halfMax-y0))/(y1-y0);
    x0 = idxAft-1;
    x1 = idxAft;
    y0 = deltaFTraces(idxAft-1,ii);
    y1 = deltaFTraces(idxAft,ii);
    idxAftInt(ii) = (x0*(y1-halfMax)+x1*(halfMax-y0))/(y1-y0);
    
    eventFWHM(ii) = (idxAftInt(ii)-idxBefInt(ii))*infoFCM.dt;
end
% Remove very large numbers
eventFWHM(eventFWHM>mean(eventFWHM*2)) = [];

hfigFWHM = figure(20); clf;
fwhmBinWidth = 0.5;
fwhmEdges = (0:fwhmBinWidth:ceil(max(eventFWHM)));
histFWHM = histogram(eventFWHM,fwhmEdges);

title({sprintf('Mean FWHM = %.2f s',mean(eventFWHM));...
    sprintf('Std FWHM = %.2f s',std(eventFWHM,1))});
histFWHM.FaceColor = [0.5 0.5 0.5]; box off;

set(hfigFWHM,'Color',[1 1 1]);
set(hfigFWHM,'units','inches','position',[9,2,1*sc,1.4*sc]);
set(gca,'tickdir','Out','ticklength',[0.02 0.02],...
    'Xtick',[0 ceil(max(eventFWHM))/2 ceil(max(eventFWHM))],...
    'YTick',[0 max(histFWHM.BinCounts)],...
    'XLim',[0 ceil(max(eventFWHM))],...
    'YLim',[0 max(histFWHM.BinCounts)]);
xlabel('Event duration (s)');
ylabel('Number of events');

set(gca,'FontSize',fontsize);
 
print([pnameOut,'/Outputs/',fnameOut,'_FWHM'],...
    '-dtiff', dpiOut);


%% ROC Curve

%Calculate roc
roc_data=[];
roc_data(1:numEvents,1)=deltaFPre;
roc_data(1:numEvents,2)=zeros(numEvents,1);
roc_data(1+numEvents:2*numEvents,1)=deltaFMax;
roc_data(1+numEvents:2*numEvents,2)=ones(numEvents,1);
roc=roc_calc(roc_data);
sc=5;

hfigROC = figure(21); clf;
set(hfigROC,'Color',[1 1 1]);
set(hfigROC,'units','inches','position',[9,1,1*sc,1.2*sc]);
plot(zscore(roc.xf),zscore(roc.yf),'-b','linewidth',1); grid on;

figure(22); clf;
hold on;
h1 = histogram(deltaFPre);
h2 = histogram(deltaFMax);
h1.BinWidth = 0.05;
h2.BinWidth = 0.05;
h1.FaceColor = 'b';
h2.FaceColor = 'r';


