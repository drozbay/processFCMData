function [ infoStruct ] = getTimeSeriesInfo( fullPathToSeries )
% GETTIMESERIESINFO
% Extracts various information from a time series tif
% stack and saves it to a structure
% Inputs:
% fullPathToSeries - String that specifies a .tif image stack of a movie or
% time series. Frame spacing must be uniform. Metadata should be saved in
% imageJ format. Necessary field is framerate (as 'finterval' or 'fps').
% Outputs:
% infoStruct - Struct of all extracted parameters

infos = imfinfo(fullPathToSeries);
% Add path to struct
infoStruct.fullPath = fullPathToSeries;
% Extract simple info from metadata
infoStruct.numFrames = numel(infos);
infoStruct.xPixels = infos(1).Width;
infoStruct.yPixels = infos(1).Height;
infoStruct.bitDepth = infos(1).BitDepth;
% Get frame interval
infoText = textscan(infos(1).ImageDescription,'%s %s','Delimiter','=');
idx = find(cellfun('length',regexp(infoText{1},'finterval'))==1);
if length(idx)==1
    finterval = str2double(infoText{2}{idx});
else
    idx = find(cellfun('length',regexp(infoText{1},'fps'))==1);
    if length(idx)==1
        finterval = 1/str2double(infoText{2}{idx});
    else
        error('Could not establish frame time');
    end
end
infoStruct.dt = finterval;
% Calculate duration
infoStruct.duration = finterval*(numel(infos)-1);
% Get pixel resolution
pixelRes = infos(1).XResolution;
if isempty(pixelRes)
    warning('No pixel size found, defaulting to 1');
    infoStruct.pixelRes = 1;
else
    infoStruct.pixelRes = pixelRes;
end

% Get resolution unit
idx = find(cellfun('length',regexp(infoText{1},'unit'))==1);
if length(idx)==1
    resUnit = infoText{2}{idx};
    if strcmpi(resUnit,'none')
        warning('Could not find resolution unit, default to mm');
        infoStruct.resUnit = 'mm';
    else
        infoStruct.resUnit = sprintf(strrep(resUnit,'\u','\x'));
    end
else
    warning('Could not find resolution unit, default to mm');
    infoStruct.resUnit = 'mm';
end

end

