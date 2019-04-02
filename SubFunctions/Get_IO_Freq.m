function [MAG,PHASE] = Get_IO_Freq(Freq,Mag,Phase,uFreq,varargin)
%% Get_IO_Freq: gets magnitude & phase values at specified frequenices
%   INPUTS:
%       Freq  	: frequency vector
%       Mag   	: magnitude vector
%       Phase 	: Phase vector
%       uFreq   : unique freuncies present
%       fTol    : frequency search tolerance (default=2*T) >>> varargin{1}
%       debug 	: showplot (default=false) >>> varargin{2}
%   OUTPUTS:
%       MAG     : manitude at uFreq
%       PHASE  	: phase at uFreq
%---------------------------------------------------------------------------------------------------------------------------------
if nargin==4
    fTol = 2*mean(diff(Freq)); % default tolerance range
    debug = false;
elseif nargin==5
    fTol = varargin{1};
 	debug = false;
elseif nargin==6
    fTol = varargin{1};
    debug = varargin{2};
elseif nargin>6
    error('Too many inputs')
elseif nargin<4
    error('Not enough inputs: requires four')
else
    error('Not Possible')
end

nFreq = length(uFreq); % # of frequencies to use
MAG = nan(nFreq,1);
PHASE = nan(nFreq,1);
for kk = 1:nFreq
    fRange = [uFreq(kk)-fTol , uFreq(kk)+fTol]; % frequency search range
    idxRange = Freq>=fRange(1) & Freq<=fRange(2); % index search range
    magRange = Mag(idxRange); % magnitude search range
    fIdx = Mag==max(magRange); % index of max magnitude
    
    MAG(kk) = max(magRange); % store magnitude at uFreq
    PHASE(kk) = Phase(fIdx); % store phase at uFreq
    
    if debug % plot magnitude & phase
    figure ; clf
        subplot(2,1,1) ; hold on ; ylabel('Magnitude')
        plot(Freq,Mag,'*-b')
        plot(uFreq(kk),MAG(kk),'or')
        xlim([0 max(uFreq+1)])

        subplot(2,1,2) ; hold on ; ylabel('Phase')
        plot(Freq,Phase,'*-b')
        plot(uFreq(kk),PHASE(kk),'or')
        xlim([0 max(uFreq+1)])
        xlabel('Time')
    end
end
end