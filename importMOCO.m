function [ mocoData ] = importMOCO( fullPath, infoStruct )
%% Initialize variables.
delimiter = ',';
startRow = 2;
endRow = inf;

% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
formatSpec = '%f%f%f%[^\n\r]';

% Open the text file.
fileID = fopen(fullPath,'r');

% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

% Close the text file.
fclose(fileID);

% Rescale data to actual dimensions
dt = infoStruct.dt;
dx = 1/infoStruct.pixelRes;
dy = dx;
timeArray = dataArray{1}*dt-dt;
xArray = dataArray{2}*dx;
yArray = dataArray{3}*dy;
magArray = sqrt(xArray.^2+yArray.^2);

%% Create output variable
mocoData = table(timeArray,xArray,yArray,magArray, 'VariableNames', {'t_s','x_um','y_um','mag_um'});


end

