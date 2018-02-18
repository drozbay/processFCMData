function [ imMaxProjActive, imFilt_wEvent_short, imFilt_noEvent_short ] = getActiveFrames( imInput, activeROIData, infoStruct, rWin )
%GETACTIVEFRAMES
%% Filter image stack
numFrames = infoStruct.numFramesTrim;
frameStart = infoStruct.trimFrames(1);
frameEnd = infoStruct.trimFrames(2);
% Trim input image
imInputTrim = imInput(:,:,frameStart:frameEnd);

imFilt = zeros(size(imInputTrim));
for ii=1:numFrames
    if rWin>1
        imFilt(:,:,ii) = imgaussfilt(medfilt3(imInputTrim(:,:,ii),[1,1,rWin]),1);
    else
        imFilt(:,:,ii) = imgaussfilt(medfilt3(imInputTrim(:,:,ii),[1,1,1]),1);
    end
end

%% Isolate frames with events
imFilt_noEvent = imFilt;
imFilt_wEvent = imFilt;
% Clear contents of frames that have events
imFilt_noEvent = padarray(imFilt_noEvent,[0 0 rWin],0);
numEvents = length(horzcat(activeROIData.eventMaxIdx));
eventIdx = sort(horzcat(activeROIData.eventMaxIdx))+rWin;
for jj=1:numEvents
    imFilt_noEvent(:,:,eventIdx(jj)-rWin:eventIdx(jj)+rWin) = ...
        zeros(size(imFilt,1),size(imFilt,2),rWin*2+1);
end
% Image stack with frames that have events cleared
imFilt_noEvent = imFilt_noEvent(:,:,rWin+1:end-rWin);
% Indices of frames that do not have events
noEventIdx = squeeze(max(max(imFilt_noEvent,[],1),[],2))>0;
% Image stack with frames that don't have events cleared
imFilt_wEvent(:,:,noEventIdx) = zeros(size(imFilt,1),size(imFilt,2),...
    sum(noEventIdx));
% Image stacks with all blank frames removed
imFilt_noEvent_short = imFilt(:,:,noEventIdx);
imFilt_wEvent_short = imFilt(:,:,~noEventIdx);

% Average image of all frames without events
imFilt_noEvent_mean = mean(imFilt_noEvent_short,3);
% Maximum projection of 
imMaxProjActive = max(imFilt_wEvent_short-imFilt_noEvent_mean,[],3);

end

