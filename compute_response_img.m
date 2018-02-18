%This code determines the pixels in a movie where there are responses
%larger than times_sdx median SD of the intensity

clear all;
close all;

w=11;  %Number of pixels for the round window used to determine response areas
twin=5; %Time window
times_sd=5; %criterion for response

min_size=10; %Area threshold in pixels

% Get image path and filename from user
[fname,pname,nCancel] = uigetfile({'*.tif;*.tiff'},'Select the TIMELAPSE file...');
if nCancel
    inputPath = [pname,fname];
    pnameStart = pname;
    %save('timeCorr_cfg.mat','pnameStart','-append');
else
    error('Cancelled')
end

%% Import images
% Get image information
infoInput = imfinfo(inputPath);
% Number of frames to drop
frameDropStart = 0;
frameDropEnd = 200;
% numImages is final frame to collect
numImages = length(infoInput)-frameDropStart-frameDropEnd;
xInput = infoInput.Width;
yInput = infoInput.Height;
% Set up matrix
imInput = double(zeros(yInput,xInput,numImages));
% Read images into stack
hWait = waitbar(0, sprintf('Importing %d images...',numImages));
for ii = 1:numImages
    imInput(:,:,ii) = double(imread(inputPath,'tif',ii+frameDropStart));
    waitbar(ii/numImages,hWait);
end
close(hWait);

%% Filter and remove average background
imInput_MedGauss = zeros(size(imInput));
hWait = waitbar(0, sprintf('Processing %d images...',numImages));
for ii=1:numImages
    imInput_MedGauss(:,:,ii) = imgaussfilt(medfilt2(imInput(:,:,ii),[1,1]),1);
    waitbar(ii/numImages,hWait);
end
close(hWait);
imInput_MedGaussAverage = mean(imInput_MedGauss,3);
imInput_SubMean = imInput_MedGauss - imInput_MedGaussAverage;

% %% Show an image
% maxVal = max(imInput_SubMean(:));
% figure(1);
% imshow(imInput_SubMean(:,:,115)/maxVal);

%% Find the responses
tic
resp=zeros(xInput,yInput);
hWait = waitbar(0, sprintf('Finding ROIs...',numImages));
for y=1+(w-1)/2:yInput-(w-1)/2
    for x=1+(w-1)/2:xInput-(w-1)/2
        % Neighborhood
        a = imInput_SubMean(y-(w-1)/2:y+(w-1)/2,x-(w-1)/2:x+(w-1)/2,:);         % Extract the neighborhood
        b=zeros(1,numImages);
        b(:,:) = mean(mean(a,1),2); % Get its mean
        
        %Get the meadian std
        these_std_vals=[];
        for ii=1:numImages-twin
            these_std_vals=[these_std_vals std(b(ii:ii+twin))];
        end
        med_std=median(these_std_vals);
        
         %Calculate responses
        for ii=1:numImages-twin
            if b(ii+twin) >times_sd*med_std+mean(b(ii:ii+twin-1)) 
                resp_t(y,x,ii+twin)=1;
                resp(y,x)=1;
            else
                resp_t(y,x,ii+twin)=0;
            end
        end
        %calculate the median baseline
        med_base(y,x)=median(b(resp_t(y,x,:)==0));
    end
    waitbar(y/(yInput-(w-1)/2),hWait);
end
close(hWait);

toc

figure(1)
imshow(resp)
title('Responsiveness')

%Save the responsiveness BW image
imwrite(resp,[pname fname(1:end-4) 'ROI.tif'])

pffft=1;
