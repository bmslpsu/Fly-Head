function [etho]=MakeEthogram(rawVid,playback)

global state

hfig = figure(1); %gen figure handles
set(hfig,'ButtonDownFcn',@clicker) %set button down function to clicker (see below)
vidData = squeeze(rawVid);
[~,~,nFrame] = size(vidData); % get dimensions  of video data
hfig.UserData=1; %use userdata as a flag.
etho=NaN(nFrame,1);
state=NaN;
tic
for n = playback:playback:nFrame
   imshow(vidData(:,:,n));
   etho((n-playback+1):n)=state;
   pause(.01)
end
toc

NaNidx = isnan(etho);
definedstate = etho(~NaNidx);
firststate = definedstate(1);
etho(NaNidx) = firststate;

imagesc(etho')
pause()

%what happens when you click on the figure
function clicker(hobject,~)
whichbutton=get(hobject, 'selectiontype'); %this gets what button was pressed
switch whichbutton % if anything really but i put this here anyways
    case 'normal' %left mouse click
        state=1; %set the object's user data to 0 (ie 0 for live feed)
        disp('walk')
    case 'alt' %right mouse click
        state=0;
        disp('stand')
end
end
end