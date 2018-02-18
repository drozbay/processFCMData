function [ bwROI, wsROI, lwROI, rgbROI] = getFCMROI( imInput, infoStruct, cfg, varargin)
%GETFCMROI
w = cfg.w;
tWin = cfg.tWin;
sdThresh = cfg.sdThresh;
minSize = cfg.minSize;
xInput = size(imInput,2);
yInput = size(imInput,1);
% If existing ROI is not included, run code to find active ROIs
switch nargin
case 3
    if isfield(infoStruct,'numFramesTrim')
        imInput = imInput(:,:,infoStruct.trimFrames(1):infoStruct.trimFrames(2));
    end
    numFrames = size(imInput,3);
    %% Filter and remove average background
    imInput_MedGauss = zeros(size(imInput));
    hWait = waitbar(0, sprintf('Removing BG from %d images...',numFrames));
    for ii=1:numFrames
        imInput_MedGauss(:,:,ii) = imgaussfilt(medfilt2(imInput(:,:,ii),[1,1]),1);
        waitbar(ii/numFrames,hWait);
    end
    close(hWait);
    imInput_MedGaussAverage = mean(imInput_MedGauss,3);
    imInput_SubMean = imInput_MedGauss - imInput_MedGaussAverage;

    %% Pad image
    wPad = (w-1)/2;
    imInput_SubMean_Padded = padarray(imInput_SubMean,[wPad,wPad],0,'both');
    xInputPad = size(imInput_SubMean_Padded,2);
    yInputPad = size(imInput_SubMean_Padded,1);

    %% Find the active regions
    activeROI=zeros(xInputPad,yInputPad,'logical');
    these_std_vals=zeros(1,numFrames-tWin);
    hWait = waitbar(0, sprintf('Finding ROIs...'));
    totIterations = length(1+(w-1)/2:yInputPad-(w-1)/2);
    tocCurrent = 0;
    jj = 0;
    for y=1+(w-1)/2:yInputPad-(w-1)/2
        jj = jj+1;
        tic;
        for x=1+(w-1)/2:xInputPad-(w-1)/2
            % Neighborhood
            a = imInput_SubMean_Padded(y-(w-1)/2:y+(w-1)/2,x-(w-1)/2:x+(w-1)/2,:);         % Extract the neighborhood
            b = squeeze(mean(mean(a,1),2)); % Get its mean

            %Get the meadian std
            for ii=1:numFrames-tWin
                these_std_vals(ii)=std(b(ii:ii+tWin));
            end
            med_std=median(these_std_vals);

             %Calculate responses
            for ii=1:numFrames-tWin
                if b(ii+tWin)>sdThresh*med_std+mean(b(ii:ii+tWin-1)) 
                    activeROI(y,x)=1;
                end
            end
        end
        tocCurrent = tocCurrent + toc;
        waitbar(jj/totIterations,hWait,...
            sprintf('Finding ROIs... %.0f sec remaining',...
            tocCurrent*(totIterations/jj-1)) );
    end
    close(hWait);
    bwROI = logical(activeROI(wPad+1:end-wPad,wPad+1:end-wPad));
case 4
    bwROI = varargin{1};
otherwise
    error('Invalid number of input arguments')
end
%% Create outputs
% ROI operations
% Watershed the distance to center map
ROIWatershed = watershed(-bwdist(~bwROI,'euclidean'),8);
% Use the watershed lines to create borders
bwROIWatershed = bwROI;
bwROIWatershed(ROIWatershed==0) = 0;
% Remove some pixels from edges to reduce influence of motion correction
rMarg = 6;
bwROIWatershed = centerPadCrop(...
    centerPadCrop(bwROIWatershed,yInput-rMarg,xInput-rMarg),...
    yInput,xInput);
% Remove all rois smaller than area threshold and fill holes inside ROIS
wsROI = imfill(bwareaopen(bwROIWatershed, minSize),'holes');
% Label ROIs
lwROI = bwlabel(wsROI);
% Make RGB ROIs
rgbROI = label2rgb(lwROI, 'hsv', 'k', 'shuffle');

%% Save ROI file
[pname, fname, ~] = fileparts(infoStruct.fullPath);
imwrite(uint8(bwROI),[pname,'/', fname, '_ROI.tif'])

end

