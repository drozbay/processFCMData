function [flow, flowMag1D] = getBehOpticalFlow(imMovie,infoStruct,subFact)
% getBehOpticalFlow
frameStart = infoStruct.trimFrames(1);
frameEnd = infoStruct.trimFrames(2);
numFrames = infoStruct.numFramesTrim;
% Trim frames
imMovieTrim = imMovie(:,:,frameStart:frameEnd);


%% Optical flow
% Sub-sample images
framerate = 1/(infoStruct.dt*subFact);
imSubSample = imMovieTrim(:,:,1:subFact:end);
numFrameSub = size(imSubSample,3);
% Use Horn-Schunck method
opticFlow = opticalFlowHS('Smoothness',2,'MaxIteration',10,'VelocityDifference',0);
% Initialize flow structure
flow = estimateFlow(opticFlow,imSubSample(:,:,1));
flow = repmat(flow,numFrameSub,1);
% Set number of frames to estimate flow between
flowWindow = 1;
% Perform flow analysis
hWait = waitbar(0, sprintf('Performing optical flow on %d images...',numFrameSub));
for ii=1:numFrameSub-flowWindow
    flow(ii) = estimateFlow(opticFlow,imSubSample(:,:,ii));
    flow(ii) = estimateFlow(opticFlow,imSubSample(:,:,ii+flowWindow));
    waitbar(ii/numFrameSub,hWait);
end
close(hWait);

imageSize_mm = infoStruct.xPixels/infoStruct.pixelRes;
mmPerPixel = imageSize_mm/size(imMovieTrim,1);
flowMag1D = zeros(numFrameSub-flowWindow,1);
for ii=1:numFrameSub-flowWindow
%     magSorted = sort(abs(flow(ii).Magnitude(:)),'descend');
%     flowMag1D(ii) = mean(magSorted(1:10))*framerate;
    flowMag1D(ii) = max(abs(flow(ii).Magnitude(:)))*framerate*mmPerPixel*100;
end
%
flowMag1DFilt = movmedian(flowMag1D,10);
flowMag1DFilt = (flowMag1DFilt-min(flowMag1DFilt));

%%
% jj=5815;
% figure(1); clf;
% subplot(1,2,1); 
% colormap(gray)
% imagesc(imSubSample(:,:,jj));hold on; 
% plot(flow(jj),'DecimationFactor',[5 5],'ScaleFactor',100);
% subplot(1,2,2);
% plot(flowMag1D); hold on;
% plot([jj jj],[0 1],'r','linewidth',1);
% title(['Current mm/s = ',num2str(flowMag1D(jj))]);

end