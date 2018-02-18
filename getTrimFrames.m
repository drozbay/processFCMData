function [infoStructOut] = getTrimFrames(imInputStack,infoStruct)

S.imageStack=makestack(imInputStack,'uint8');
S.infoStruct=infoStruct;

%% Creating the UI
% Get user's screen size
SCRSZ=get(0,'ScreenSize');                                                  
figheight=(SCRSZ(4)-100)/2;
figwidth=(SCRSZ(4)*0.9)/2;
% Inside padding
pad=10;
% Step sizes
smallstep=1/(size(S.imageStack,3)-1);
largestep=smallstep*10;

% Figure handle
S.hfig = figure('units','pixels',...                                          
    'position',[SCRSZ(3)/4 SCRSZ(4)/4 figwidth figheight],...
    'menubar','figure',...
    'name','getTrimFrames',...
    'numbertitle','off',...
    'resize','off');

% Color limits for the plotting and colorbar. Needs to be done before callbacks are assigned
S.clims=[0 max(S.imageStack(:))];

% Create UI elements
S.mainAxes = axes('units','pixels',...                                            
    'position',[3*pad 9*pad 43*pad figheight-8*pad],...
    'fontsize',10,...
    'nextplot','replacechildren',...
    'ytick',[],'xtick',[]);
% Slider for selecting current frame
S.theSlider = uicontrol('style','slide',...                                        
    'unit','pix',...                           
    'position',[2*pad pad 42*pad 2*pad],...
    'min',1,'max',size(S.imageStack,3),'val',1,...
    'SliderStep', [smallstep largestep]);
% "Current frame:" text
S.txtCurFrame = uicontrol('style','text',...                                       
    'unit','pix',...
    'position',[14*pad 3.3*pad 10*pad 2*pad],...
    'fontsize',10,...
    'string','Current frame:');
% Edit box for selecting current frame
S.edCurFrame = uicontrol('style','edit',...                                         
    'unit','pix',...
    'position',[24*pad 3.5*pad 5*pad 2*pad],...
    'fontsize',10,...
    'string','1');
% Edit box to input start of range
S.edStart = uicontrol('style','edit',...                                         
    'unit','pix',...
    'position',[4*pad 3.5*pad 6*pad 2*pad],...
    'fontsize',10,...
    'string','1');
% Edit to input end of range
S.edEnd = uicontrol('style','edit',...                                         
    'unit','pix',...
    'position',[35*pad 3.5*pad 6*pad 2*pad],...
    'fontsize',10,...
    'string',num2str(size(S.imageStack,3)));
% Colormap menu
S.cmapList={'Gray','Jet','Hot','Copper'};
S.popCmap = uicontrol('style','popupmenu',...                  
    'unit','pix',...
    'position',[40*pad figheight-5*pad 6*pad 2*pad],...
    'String', S.cmapList);
% Button to choose current frame as start of range
S.selStart = uicontrol('style','pushbutton',...
    'unit','pix',...
    'position',[2*pad 6*pad 13*pad 3*pad],...
    'fontsize',10,...
    'string','Select Start Frame');
% Button to choose current frame as end of range
S.selEnd = uicontrol('style','pushbutton',...
    'unit','pix',...
    'position',[31*pad 6*pad 13*pad 3*pad],...
    'fontsize',10,...
    'string','Select End Frame');
% Button to submit range to main program
S.butSubmitRange = uicontrol('style','pushbutton',...                         
    'unit','pix',...
    'position',[18*pad 6*pad 10*pad 3*pad],...
    'fontsize',10,...
    'string','Submit Range');

% Setting callbacks
set([S.theSlider,S.edCurFrame],'Call',{@switchframe,S});
set(S.popCmap,'Callback',{@setcmap,S});
set(S.selStart,'Callback', {@selectStartFunc,S});
set(S.selEnd,'Callback',{@selectEndFunc,S});
set([S.edStart,S.edEnd],'Call',{@roundMe,S});
set(S.butSubmitRange,'Callback',{@submitRangeFunc,S});

% Display the first frame
imagesc(squeeze(S.imageStack(:,:,1)),S.clims);
axis equal tight
%Set colormap
setcmap(S.popCmap,[],S);                                                       
colorbar;

% Wait until mainAxes is deleted, either through its callback or by
% closing the figure
waitfor(S.mainAxes);

% Check that figure wasn't prematurely closed
if ishandle(S.hfig)
    % Send new infoStruct to output
    F = get(S.butSubmitRange,'callback');
    infoStructOut = F{2}.infoStruct;
    close(S.hfig);
else
    % Default values if figure was closed
    S.infoStruct.trimFrames = [1 size(S.imageStack,3)];
    infoStructOut = S.infoStruct;
end

end

%% Callback functions
% Colormap callback
function setcmap(varargin)
[h,S] = varargin{[1,3]};
%Create a colormap cmap with 256 colors and the chosen colormap
eval(['cmap=colormap(', lower(S.cmapList{get(h,'value')}),...
    '(', num2str(S.clims(2)), '));']);
%If "jet" is chosen, set 0 to black and 255 to white
if get(h,'value')==1
    cmap(1,:)=[0 0 0];
    cmap(end,:)=[1 1 1];
%If "gray" is chosen, set 0 to blue and 255 to red
elseif get(h,'value')==2
    cmap(1,:)=[0 0 1];
    cmap(end,:)=[1 0 0];
end
colormap(S.mainAxes,cmap);
end

%% Slider & frame edit box callback
function [] = switchframe(varargin)
[h,S] = varargin{[1,3]};
% Who called?
switch h
    % Editbox:
    case S.edCurFrame
        % Get the slider's info
        sliderstate =  get(S.theSlider,{'min','max','value'});
        % The new frame number
        enteredvalue = str2double(get(h,'string'));
        
        if enteredvalue >= sliderstate{1} && enteredvalue <= sliderstate{2} %Check if the new frame number actually exists
            slidervalue=round(enteredvalue);
            set(S.theSlider,'value',slidervalue);                                   %If it does, move the slider there
        else
            set(h,'string',sliderstate{3});                                  %User tried to set slider out of range, keep value
            return
        end
    case S.theSlider                                                              % The slider called...
        slidervalue=round(get(h,'value'));                                  % Get the new slider value
        set(S.edCurFrame,'string',slidervalue)                                      % Set editbox to current slider value
end
% Show frame
imagesc(squeeze(S.imageStack(:,:,slidervalue)),S.clims);
setcmap(S.popCmap,[],S);
end
%% Callback to round off edit boxes
function []=roundMe(varargin)
h=varargin{1};
set(h,'string',round(str2double(get(h,'string'))));
end
%% Callbacks for selecting starting and ending frames
function []=selectStartFunc(varargin)
[~,S] = varargin{[1,3]};
currentFrame = round(str2double(get(S.edCurFrame,'string')));
endFrame = str2double(get(S.edEnd,'string'));
if currentFrame>=endFrame
    set(S.edStart,'string',num2str(endFrame-1));
else
    set(S.edStart,'string',num2str(currentFrame));
end

end
function []=selectEndFunc(varargin)
[~,S] = varargin{[1,3]};
currentFrame = round(str2double(get(S.edCurFrame,'string')));
startFrame = str2double(get(S.edStart,'string'));
if currentFrame<=startFrame
    set(S.edEnd,'string',num2str(startFrame+1));
else
    set(S.edEnd,'string',num2str(currentFrame));
end

end
%% Callback for submit button
function []=submitRangeFunc(varargin)
S=varargin{3};
trimStart = str2double(get(S.edStart,'string'));
trimEnd = str2double(get(S.edEnd,'string'));
S.infoStruct.trimFrames = [trimStart trimEnd];
% Save new S struct into callback
set(S.butSubmitRange,'callback',{@submitRangeFunc,S});
% Delete button to trigger resuming main function execution
delete(S.mainAxes);
end

%% Function that creates the 3D-array used
function [outstack]=makestack(arrayorfile,dataclass)                        %First input is either a filename or a 3D array, second is the class

if ischar(arrayorfile)                                                      %If it's a filename, read the file using tiffread
    lsm= tiffread(arrayorfile);
    
    if length(lsm(1,1).data)>1                                              %If there are several channels, pick which one to display
        chstrgs = inputdlg({'Pick channel to display:'},...
            'Enter channel',1,{'1'});
        ch=str2double(chstrgs{1});
    end
    outstack=zeros(size(lsm(1).data{1},1), size(lsm(1).data{1},2),...       %Preallocate memory for the image stack
        length(lsm));
    for k=1:length(lsm)
        outstack(:,:,k)=lsm(k).data{ch};                                    %Build the stack. There should be a better way to do this...
    end
else
    outstack=arrayorfile;                                                   %If the input was an array, that array is our image stack
    clear arrayorfile
end

maxval=max(outstack(:));                                                    %The maximum value of the stack for colormap scaling
currentclass=class(outstack);                                               %The data class of the data

switch dataclass                                                            %Switch on the desired data class 
    case 'uint8'                                                            %User wants uint8 (quickest and least memory needed)
        if ~strcmp(currentclass,'uint8')                                    %If it's already uint8, do nothing
            outstack=uint8(255*double(outstack)/double(maxval));            %Normalize and cast to uint8
        end
        
    case 'uint16'                                                           %User wants uint16
        if strcmp(currentclass,'uint8')                                     %If the class is uint8, this is stupid.
            disp('Displaying an 8-bit image as a 16-bit image is stupid. uint8 used.')
        elseif ~strcmp(currentclass,'uint16')                               %If it's already uint16, do nothing
            outstack=uint16(65535*double(outstack)/double(maxval));         %Normalize and cast to uint16
        end
end

end



