function [strCell] = num2strcell(Array)
%% num2strcell: converts numeric rray to cell of strings
%
%   INPUTS:
%       ARRAY   : numeric array
%
%   OUTPUTS:
%       CELL   	: cell array containing the value sin "ARRAY" converted to strings
%

cellArray   = num2cell(Array);
strCell     = cellfun(@(x) num2str(x), cellArray, 'UniformOutput', false);
end