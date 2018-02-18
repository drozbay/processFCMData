% StackSlider(varargin): GUI for displaying images from, e.g., an lsm stack
% or other 3D array. Can be called in several ways:
%
% StackSlider() without input lets the user choose an lsm file to be
% displayed. Images will be converted to and shown on an 8-bit scale.
%
% StackSlider(I) where I is a 3D array will display the array, assuming that
% consecutive images are arrayed along the third dimension. As above, uses
% 8-bit colors.
%
% StackSlider('filename') where filename is a path and name of a valid image
% stack will display said stack. As above, uses 8-bit colors.
%
% StackSlider(X,'type') where X is _either_ an array or filename as above
% will convert the input array to the class 'type' (uint8 or uint16) and
% display it. Note: StackSlider will not convert an 8-bit image to a 16-bit,
% since this is plain stupid.
%
% The GUI itself lets the user scroll back and forth among the frames using a
% slider, och jump to a frame by writing the frame number in an editbox.
% The colormap can be chosen from a popup menu.
% 
% Smoothing can be chosen as "no smoothing" (displays the raw data), "Gaussian"
% (user can set radius of filter and standard deviation of filter) and 
% "averaging" (moving average over a disk of user-controlled radius).
% NOTE: The size of the gaussian filter is the radius, and not the diameter 
% set when using, e.g., fspecial and imfilter.
%
% The "reset all" button does exactly that.
%
% The "make figure" button makes a new figure of the currently shown frame,
% using the selected smoothing options. Note that the new figure does not
% inherit the colormap limits, but is streched over the colormap using "imagesc".
%
% This program calls tiffread.m, which can be downloaded for free from
% http://www.cytosim.org/other/
%
% Written by Otto Manneberg, SciLifeLab 2011-07-22.
% otto.manneberg@scilifelab.se
% Copyright (c) 2011, Otto Manneberg
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of the Science for Life Laboratory (SciLifeLab) nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



function [infoStructOut] = stackSliderBNO(inputImage)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input argument checking and handling                                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ~= 1
    disp('Not enough input arguments');                 %Check for wrong number of inputs
    return
end

S.I=makestack(inputImage,'uint8');
S.infoStruct = infoStruct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Build the figure for the GUI.                                          
% All handles and the image stack are stored in the struct SS             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SCRSZ=get(0,'ScreenSize');                                                  %Get user's screen size
figheight=(SCRSZ(4)-100)/2;                                                     %A reasonable height for the GUI
figwidth=(SCRSZ(4)*0.9)/2;                                                      %A reasonable width for the GUI (the height of the screen*1.1)
pad=10;                                                                     %Inside padding in the GUI
smallstep=1/(size(S.I,3)-1);                                                %Step the slider will take when moved using the arrow buttons: 1 frame
largestep=smallstep*10;                                                     %Step the slider will take when moved by clicking in the slider: 10 frames

%%%%%%Create the figure itself. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S.fh = figure('units','pixels',...                                          
    'position',[SCRSZ(3)/4 SCRSZ(4)/4 figwidth figheight],...
    'menubar','figure',...
    'name','StackSlider',...
    'numbertitle','off',...
    'resize','off');

%%%%%%Create the axes for image display. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S.ax = axes('units','pixels',...                                            
    'position',[3*pad 9*pad 43*pad figheight-8*pad],...
    'fontsize',10,...
    'nextplot','replacechildren',...
    'ytick',[],'xtick',[]);
%%%%%%Create a slider and an editbox for picking frames. %%%%%%%%%%%%%%%%%%
S.sl = uicontrol('style','slide',...                                        
    'unit','pix',...                           
    'position',[2*pad pad 42*pad 2*pad],...
    'min',1,'max',size(S.I,3),'val',1,...
    'SliderStep', [smallstep largestep]);
S.cmtext=uicontrol('style','text',...                                       
    'unit','pix',...
    'position',[14*pad 3.3*pad 10*pad 2*pad],...
    'fontsize',10,...
    'string','Current frame:');
S.ed = uicontrol('style','edit',...                                         
    'unit','pix',...
    'position',[24*pad 3.5*pad 5*pad 2*pad],...
    'fontsize',10,...
    'string','1');
S.editStart = uicontrol('style','edit',...                                         
    'unit','pix',...
    'position',[4*pad 3.5*pad 6*pad 2*pad],...
    'fontsize',10,...
    'string','1');
S.editEnd = uicontrol('style','edit',...                                         
    'unit','pix',...
    'position',[35*pad 3.5*pad 6*pad 2*pad],...
    'fontsize',10,...
    'string',num2str(size(S.I,3)));

% Colormap menu
S.cmstr={'Gray','Jet','Hot','Copper'};                                      %Strings with the allowed colormaps
S.cmpopup = uicontrol('style','popupmenu',...                               %Popup menu for picking                        
    'unit','pix',...
    'position',[40*pad figheight-5*pad 6*pad 2*pad],...
    'String', S.cmstr);
% Button to choose current frame as start of range
S.selectStart = uicontrol('style','pushbutton',...
    'unit','pix',...
    'position',[2*pad 6*pad 13*pad 3*pad],...
    'fontsize',10,...
    'string','Select Start Frame');
% Button to choose current frame as end of range
S.selectEnd = uicontrol('style','pushbutton',...
    'unit','pix',...
    'position',[31*pad 6*pad 13*pad 3*pad],...
    'fontsize',10,...
    'string','Select End Frame');
% Button to submit range to main program
S.submitRange = uicontrol('style','pushbutton',...                         
    'unit','pix',...
    'position',[18*pad 6*pad 10*pad 3*pad],...
    'fontsize',10,...
    'string','Submit Range');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw the first frame of the stack and set callback functions           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%Draw the first frame of the stack%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S.Clims=[0 max(S.I(:))];                                                    %Color limits for the plotting and colorbar. Needs to be done before callbacks are assigned
imagesc(squeeze(S.I(:,:,1)),S.Clims);                                       %Display the first frame
axis equal tight                                                            %Make sure it's to scale
setcm(S.cmpopup,[],S)                                                       %Set colormap
colorbar                                                                    %Display a colorbar

%%%%%%Set callback functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set([S.ed,S.sl],'call',{@switchframe,S});                                   %Shared callback function for fram selection slider and editbar
set(S.cmpopup,'Callback',{@setcm,S});                                       %Callback function for changing colormap
set(S.selectStart,'Callback', {@selectStartFunc,S});                        %Callback function for the reset button
set(S.selectEnd,'Callback',{@selectEndFunc,S});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Change colormap callback function                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=setcm(varargin)                                                 %varargin is {calling handle, eventdata, struct SS}, where eventdata is empty (currently unused) when called as callback
[h,S] = varargin{[1,3]};                                                            %Extract handle of calling object
eval(['cmap=colormap(', lower(S.cmstr{get(h,'value')}),...                 %Create a colormap cmap with 256 colors and the chosen colormap
    '(', num2str(S.Clims(2)), '));']);
if get(h,'value')==1                                                        %If "jet" is chosen, set 0 to black and 255 to white
    cmap(1,:)=[0 0 0];
    cmap(end,:)=[1 1 1];
elseif get(h,'value')==2                                                    %If "gray" is chosen, set 0 to blue and 255 to red (like the range indicator on a Zeiss mic)
    cmap(1,:)=[0 0 1];
    cmap(end,:)=[1 0 0];
end
    
colormap(S.ax,cmap);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Move slider or write in frame editbox callback function                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = switchframe(varargin)                                         %varargin is {calling handle, eventdata, struct S}, where eventdata is empty (currently unused) when called as callback
[h,S] = varargin{[1,3]};                                                    %Extract handle of calling object and the struct S

switch h                                                                    %Who called?
    case S.ed                                                               %The editbox called...
        sliderstate =  get(S.sl,{'min','max','value'});                     % Get the slider's info
        enteredvalue = str2double(get(h,'string'));                         % The new frame number
        
        if enteredvalue >= sliderstate{1} && enteredvalue <= sliderstate{2} %Check if the new frame number actually exists
            slidervalue=round(enteredvalue);
            set(S.sl,'value',slidervalue);                                   %If it does, move the slider there
        else
            set(h,'string',sliderstate{3});                                  %User tried to set slider out of range, keep value
            return
        end
    case S.sl                                                               % The slider called...
        slidervalue=round(get(h,'value'));                                  % Get the new slider value
        set(S.ed,'string',slidervalue)                                      % Set editbox to current slider value
end

imagesc(squeeze(S.I(:,:,slidervalue)),S.Clims);
setcm(S.cmpopup,[],S);

end

function []=selectStartFunc(varargin)                                       % varargin is {calling handle, eventdata, struct S}, where eventdata is empty (currently unused) when called as callback
S=varargin{3};                                                              % Extract the struct S

end

function []=selectEndFunc(varargin)                                        % varargin is {calling handle, eventdata, struct S}, where eventdata is empty (currently unused) when called as callback
S=varargin{3};                                                             % Extract the struct S

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



