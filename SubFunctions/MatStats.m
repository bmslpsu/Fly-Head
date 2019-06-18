function [Stats] = MatStats(matData,dim)
%% MatStats: statistics of matrix
%   INPUTS:
%       matData  	:   matrix data
%       dim         :   dimension
%   OUTPUTS:
%       Stats       :   statistics object
%---------------------------------------------------------------------------------------------------------------------------------
Stats.Mean      = nanmean(matData,dim);
Stats.Median    = nanmedian(matData,dim);
Stats.STD       = nanstd(matData,[],dim);

% if dim==1
%     Stats.First = matData(:,1);
%     Stats.Last  = matData(:,end);
% elseif dim==2
% 	Stats.First = matData(1,:);
%    	Stats.Last  = matData(end,:);
% end


end