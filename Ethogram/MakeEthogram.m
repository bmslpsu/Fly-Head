
function [etho]=MakeEthogram(rawVid,playback)
% Allows user to manually track the flys state as walking or standing. 
% Click in the figure margins. Left click when the fly starts walking, 
% right click when the fly starts standing. Outputs a row vector with 
% 1 (walking) and 0 (stamding) for each frame of the video. 

global state                        % sets state as global to be used in subfunction

hfig = figure(1); % figure handles
set(hfig,'ButtonDownFcn',@clicker)  % mouse click calls clicker subfunction
vidData = squeeze(rawVid);          % removes extraneous video dimension
[~,~,nFrame] = size(vidData);       % get # of frames in video 
hfig.UserData=1;                    % use userdata as a flag.
etho=NaN(nFrame,1);                 % preallocates etho as NaN vector
state=NaN;                          % state is initally NaN
tic
for n = playback:playback:nFrame    % loops through frames in playback increments
   imshow(vidData(:,:,n));
   etho((n-playback+1):n)=state;    % sets current value and previous to state
   pause(.01)
end
toc

nanIdx = isnan(etho);               % sets all initial values to first declared state
stateIdx = etho(~nanIdx);
etho(nanIdx)=stateIdx(1);


imagesc(etho')                      % shows results of tracking
pause(1)


function clicker(hobject,~)         % called on button click
whichbutton=get(hobject, 'selectiontype'); % which button was pressed
switch whichbutton 
    case 'normal'                   % left mouse click
        state=1;                    % sets state to 1
        disp('walk')
    case 'alt'                      % right mouse click
        state=0;                    % sets state to 0
        disp('stand')
end
end
end