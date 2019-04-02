function [COHR] = Get_IO_Cohr(Freq,Mag,uFreq,varargin)
%% Get_IO_Cohr: gets coherence values at specified frequenices
%   INPUTS:
%       Freq  	: frequency vector
%       Mag   	: magnitude vector
%       uFreq   : unique freuncies present
%       fTol    : frequency search tolerance (default=2*T) >>> varargin{1}
%       debug 	: showplot (default=false) >>> varargin{2}
%   OUTPUTS:
%       COHR     : coherence at uFreq
%---------------------------------------------------------------------------------------------------------------------------------
if nargin==3
    fTol = 2*mean(diff(Freq)); % default tolerance range
    debug = false;
elseif nargin==4
    fTol = varargin{1};
 	debug = false;
elseif nargin==5
    fTol = varargin{1};
    debug = varargin{2};
elseif nargin>5
    error('Too many inputs')
elseif nargin<3
    error('Not enough inputs: requires four')
else
    error('Not Possible')
end

nFreq = length(uFreq); % # of frequencies to use
COHR = nan(nFreq,1);
for kk = 1:nFreq
    fRange = [uFreq(kk)-fTol , uFreq(kk)+fTol]; % frequency search range
    idxRange = Freq>=fRange(1) & Freq<=fRange(2); % index search range
    magRange = Mag(idxRange); % magnitude search range
%     fIdx = Mag==max(magRange); % index of max magnitude
    
    COHR(kk) = max(magRange); % store magnitude at uFreq
    
    if debug % plot coherence
    figure ; clf ; hold on ; ylabel('Magnitude')
        plot(Freq,Mag,'*-b')
        plot(uFreq(kk),COHR(kk),'or')
        xlim([0 max(uFreq+1)])
    end
end
end