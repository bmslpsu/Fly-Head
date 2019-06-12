function [MAG,PHASE,FIdx] = Get_IO_Freq(Freq,Mag,Phase,uFreq,varargin)
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
%       fIdx    : index at UFreq
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
nn      = length(Freq);
nFreq   = length(uFreq);    % # of frequencies to use
Ff      = nan(nFreq,1);     % closest frequency that max mag occurs
MAG     = nan(nFreq,1);     % magnitude at frequencies
PHASE   = nan(nFreq,1);     % phase at frequencies
FIdx    = nan(nn,nFreq);    % indicies at frequencies
for kk = 1:nFreq
    fRange = [uFreq(kk)-fTol , uFreq(kk)+fTol]; % frequency search range
    idxRange = Freq>=fRange(1) & Freq<=fRange(2); % index search range
    magRange = Mag(idxRange); % magnitude search range
    fIdx = Mag==max(magRange); % index of max magnitude
    Ff(kk) = Freq(fIdx);
    
    MAG(kk) = max(magRange); % store magnitude at uFreq
    PHASE(kk) = Phase(fIdx); % store phase at uFreq
    
    FIdx(:,kk) = fIdx;
end

FIdx = logical(sum(FIdx,2));
if sum(FIdx,1)<nFreq
    error('Error: conflicting frequencies')
end
[FIdx,~] = find(FIdx==true);


if debug % plot magnitude & phase
figure ; clf
    subplot(2,1,1) ; hold on ; ylabel('Magnitude')
    plot(Freq,Mag,'*-b')
    plot(Ff,MAG,'or')
    xlim([0 max(uFreq+1)])

    subplot(2,1,2) ; hold on ; ylabel('Phase')
    plot(Freq,Phase,'*-b')
    plot(Ff,PHASE,'or')
    xlim([0 max(uFreq+1)])
    xlabel('Time')
end
end